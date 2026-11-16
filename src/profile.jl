# Profiles and profile–profile alignment, ported from MUSCLE 3.8.31:
# ProfileFromMSA (profilefrommsa.cpp, msa2.cpp), the nucleotide
# sum-of-pairs column score ScoreProfPos2SPN (scorepp.cpp) and the
# dynamic programming NWSmall with its traceback (nwsmall.cpp). Scores are
# Float32 and are summed in MUSCLE's order, so paths and ties match.

# MUSCLE's nucleotide matrix NUC_SP (nucmx.cpp): BLASTZ scores
# (Chiaromonte et al. 2002) plus a centre of 60.
const _NUC_SP = Float32[91 -114 -31 -123; -114 100 -125 -31; -31 -125 100 -114; -123 -31 -114 91] .+ 60.0f0
const _GAPOPEN = -400.0f0            # SetDefaultsSPN_DNA
const _MINUS_INF = -1f37

struct _Profile
    counts::Matrix{Float32}       # 4 × L weighted letter counts
    order::Matrix{Int8}           # 4 × L letters by decreasing count
    scores::Matrix{Float32}       # 4 × L: Σ_j counts[j]·NUC_SP[i, j]
    gapopen::Vector{Float32}
    gapclose::Vector{Float32}
end

# rows of `msa` are ASCII letters ('-' gaps); `w` their weights (sum 1)
function _profile(msa::AbstractMatrix{UInt8}, w::AbstractVector{Float32}, sc)
    m, L = size(msa)
    P = _Profile(zeros(Float32, 4, L), zeros(Int8, 4, L), zeros(Float32, 4, L),
                 zeros(Float32, L), zeros(Float32, L))
    for c in 1:L
        cnt = view(P.counts, :, c)
        for r in 1:m
            x = msa[r, c]
            x == _GAP && continue
            _addcounts!(cnt, _letter(x), w[r])
        end
        _sortorder!(view(P.order, :, c), cnt)
        for i in 1:4
            s = 0.0f0
            for j in 1:4
                s += cnt[j] * sc.matrix[i, j]
            end
            P.scores[i, c] = s
        end
        start = stop = 0.0f0
        for r in 1:m
            msa[r, c] == _GAP || continue
            (c == 1 || msa[r, c-1] != _GAP) && (start += w[r])
            (c == L || msa[r, c+1] != _GAP) && (stop += w[r])
        end
        P.gapopen[c] = Float32(1.0 - start) * sc.gapopen / 2
        P.gapclose[c] = Float32(1.0 - stop) * sc.gapopen / 2
    end
    return P
end

# a letter's weight into the counts (GetFractionalWeightedCounts). For
# nucleotides MUSCLE's switch tests amino-acid codes, so only X (whose
# nucleotide code equals the amino-acid code of R) is split between G and
# A; every other wildcard adds w/20 to each base.
@inline function _addcounts!(cnt, l, w)
    if l < 4
        cnt[l+1] += w
    elseif l == _NX_X
        cnt[_NX_G+1] += w / 2
        cnt[_NX_A+1] += w / 2
    else
        f = w / 20
        for b in 1:4
            cnt[b] += f
        end
    end
end

# bubble sort of the letters by decreasing count (SortCounts)
function _sortorder!(order, cnt)
    for i in 1:4
        order[i] = i
    end
    swapped = true
    while swapped
        swapped = false
        for k in 1:3
            a, b = order[k], order[k+1]
            if cnt[a] < cnt[b]
                order[k], order[k+1] = b, a
                swapped = true
            end
        end
    end
end

# ScoreProfPos2SPN: A's letters in decreasing order of count, until a zero
@inline function _colscore(A::_Profile, i, B::_Profile, j)
    s = 0.0f0
    @inbounds for k in 1:4
        l = A.order[k, i]
        f = A.counts[l, i]
        f == 0 && break
        s += f * B.scores[l, j]
    end
    return s
end

# TERMGAPS_Half: no penalty for opening a gap before the first column or
# closing one after the last
function _termgaps!(P::_Profile)
    L = length(P.gapopen)
    L == 0 && return
    P.gapopen[1] != _MINUS_INF && (P.gapopen[1] = 0.0f0)
    L > 1 && P.gapclose[L] != _MINUS_INF && (P.gapclose[L] = 0.0f0)
end

const _BIT_DM, _BIT_IM, _BIT_xM, _BIT_MD, _BIT_MI = 0x01, 0x02, 0x03, 0x04, 0x08

# NWSmall: the best global path through profiles A and B, as a vector of
# 'M' (column of A with column of B), 'D' (column of A against gaps) and
# 'I' (column of B against gaps).
function _nwsmall(A::_Profile, B::_Profile, e::Float32)
    LA, LB = length(A.gapopen), length(B.gapopen)
    LA == 0 && return fill(UInt8('I'), LB)
    LB == 0 && return fill(UInt8('D'), LA)
    _termgaps!(A)
    _termgaps!(B)
    NI = _MINUS_INF
    TB = zeros(UInt8, LA + 1, LB + 1)          # TB[i+1, j+1]: cell (i, j)
    MPrev = fill(NI, LB + 1); MCurr = fill(NI, LB + 1); MNext = fill(NI, LB + 1)
    DRow = fill(NI, LB + 1)                    # index j+1: column j
    MPrev[1] = 0.0f0
    MCurr[1] = NI
    MCurr[2] = _colscore(A, 1, B, 1)
    for j in 2:LB
        MCurr[j+1] = _colscore(A, 1, B, j) + B.gapopen[1] + (j - 2) * e + B.gapclose[j-1]
        TB[2, j+1] = (TB[2, j+1] & ~_BIT_xM) | _BIT_IM
    end
    Iij = NI
    for i in 1:LA-1
        Iij = NI
        DRow[1] = A.gapopen[1] + (i - 1) * e
        MCurr[1] = NI
        if i == 1
            MCurr[2] = _colscore(A, 1, B, 1)
            TB[2, 2] &= ~_BIT_xM
        else
            MCurr[2] = _colscore(A, i, B, 1) + A.gapopen[1] + (i - 2) * e + A.gapclose[i-1]
            TB[i+1, 2] = (TB[i+1, 2] & ~_BIT_xM) | _BIT_DM
        end
        for j in 1:LB-1
            MNext[j+2] = _colscore(A, i + 1, B, j + 1)
        end
        for j in 1:LB-1
            Iij = _recurse_d_i!(TB, DRow, MPrev, MCurr, Iij, A, B, i, j, e)
            # M at (i+1, j+1)
            DM = DRow[j+1] + A.gapclose[i]
            IM = Iij + B.gapclose[j]
            MM = MCurr[j+1]
            t = TB[i+2, j+2] & ~_BIT_xM
            if MM >= DM && MM >= IM
                MNext[j+2] += MM
            elseif DM >= MM && DM >= IM
                MNext[j+2] += DM
                t |= _BIT_DM
            else
                MNext[j+2] += IM
                t |= _BIT_IM
            end
            TB[i+2, j+2] = t
        end
        Iij = _recurse_d_i!(TB, DRow, MPrev, MCurr, Iij, A, B, i, LB, e)
        MPrev, MCurr, MNext = MCurr, MNext, MPrev
    end
    MCurr[1] = NI
    if LA > 1
        MCurr[2] = _colscore(A, LA, B, 1) + (LA - 2) * e + A.gapopen[1] + A.gapclose[LA-1]
    else
        MCurr[2] = _colscore(A, LA, B, 1) + A.gapopen[1] + A.gapclose[1]
    end
    TB[LA+1, 2] = (TB[LA+1, 2] & ~_BIT_xM) | _BIT_DM
    DRow[1] = NI
    for j in 1:LB
        _recurse_d!(TB, DRow, MPrev, A, LA, j, e)
    end
    Iij = NI
    for j in 1:LB
        Iij = _recurse_i!(TB, MCurr, Iij, B, LA, j, e)
    end
    score, edge = MCurr[LB+1], UInt8('M')
    DRow[LB+1] > score && ((score, edge) = (DRow[LB+1], UInt8('D')))
    Iij > score && ((score, edge) = (Iij, UInt8('I')))
    return _traceback(TB, LA, LB, edge)
end

@inline function _recurse_d_i!(TB, DRow, MPrev, MCurr, Iij, A, B, i, j, e)
    _recurse_d!(TB, DRow, MPrev, A, i, j, e)
    return _recurse_i!(TB, MCurr, Iij, B, i, j, e)
end

# D(i, j) = max(D(i−1, j) + e, M(i−1, j) + gap open of A's column i)
@inline function _recurse_d!(TB, DRow, MPrev, A, i, j, e)
    DD = DRow[j+1] + e
    MD = MPrev[j+1] + A.gapopen[i]
    if DD > MD
        DRow[j+1] = DD
    else
        DRow[j+1] = MD
        TB[i+1, j+1] |= _BIT_MD
    end
end

# I(i, j) = max(I(i, j−1) + e, M(i, j−1) + gap open of B's column j)
@inline function _recurse_i!(TB, MCurr, Iij, B, i, j, e)
    Iij += e
    MI = MCurr[j] + B.gapopen[j]
    if MI >= Iij
        Iij = MI
        TB[i+1, j+1] |= _BIT_MI
    end
    return Iij
end

# BitTraceBack
function _traceback(TB, LA, LB, last::UInt8)
    path = UInt8[]
    sizehint!(path, LA + LB)
    i, j = LA, LB
    edge = last
    while true
        push!(path, edge)
        bits = TB[i+1, j+1]
        next = edge == UInt8('M') ? (bits & _BIT_xM == _BIT_DM ? UInt8('D') : bits & _BIT_xM == _BIT_IM ? UInt8('I') : UInt8('M')) :
               edge == UInt8('D') ? (bits & _BIT_MD != 0 ? UInt8('M') : UInt8('D')) :
                                    (bits & _BIT_MI != 0 ? UInt8('M') : UInt8('I'))
        edge == UInt8('M') ? (i -= 1; j -= 1) : edge == UInt8('D') ? (i -= 1) : (j -= 1)
        (i == 0 && j == 0) && break
        edge = next
    end
    return reverse!(path)
end

# the merged alignment of `a` and `b` along `path` (rows of `a` first)
function _mergemsa(a::Matrix{UInt8}, b::Matrix{UInt8}, path::Vector{UInt8})
    ma, mb = size(a, 1), size(b, 1)
    out = Matrix{UInt8}(undef, ma + mb, length(path))
    i = j = 0
    for (c, op) in enumerate(path)
        if op == UInt8('I')
            out[1:ma, c] .= _GAP
        else
            i += 1
            out[1:ma, c] .= view(a, :, i)
        end
        if op == UInt8('D')
            out[ma+1:end, c] .= _GAP
        else
            j += 1
            out[ma+1:end, c] .= view(b, :, j)
        end
    end
    return out
end

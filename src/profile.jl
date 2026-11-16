# Profile–profile alignment (Gotoh, affine gaps) for the progressive merges.
#
# A group of aligned sequences is a matrix of letters (sequences × columns).
# Its profile is the mean A/C/G/T content of each column plus the fraction
# of sequences with a letter there (occupancy). Two columns score the mean
# substitution score over all pairs of letters; a column set against a gap
# pays the gap penalties times its occupancy, so columns that are mostly
# gaps already are cheap to gap again. Gaps at either end pay no opening
# penalty.

struct _Scoring
    S::Matrix{Float32}     # 4 × 4 substitution scores, A C G T
    open::Float32          # gap opening penalty (> 0)
    extend::Float32        # gap extension penalty, per column (> 0)
end

function _Scoring(match::Real, mismatch::Real, gapopen::Real, gapextend::Real)
    gapopen >= 0 && gapextend >= 0 || throw(ArgumentError("gap penalties must be ≥ 0"))
    S = fill(Float32(mismatch), 4, 4)
    for a in 1:4
        S[a, a] = match
    end
    return _Scoring(S, gapopen, gapextend)
end

struct _Profile
    freq::Matrix{Float32}  # 4 × columns
    occ::Vector{Float32}
end

function _Profile(cols::Matrix{UInt8})
    m, L = size(cols)
    freq = zeros(Float32, 4, L)
    occ = zeros(Float32, L)
    for c in 1:L, r in 1:m
        x = cols[r, c]
        x == _GAP && continue
        occ[c] += 1
        for a in 1:4
            freq[a, c] += _FRAC[a, x+1]
        end
    end
    freq ./= m
    occ ./= m
    return _Profile(freq, occ)
end

const _M, _X, _Y = 0x00, 0x01, 0x02   # match; column of A vs gap; column of B vs gap

# Align the groups `a` and `b` (letter matrices); returns the merged group,
# rows of `a` first.
function _alignprofiles(a::Matrix{UInt8}, b::Matrix{UInt8}, sc::_Scoring)
    path = _gotoh(_Profile(a), _Profile(b), sc)
    ma, mb = size(a, 1), size(b, 1)
    out = Matrix{UInt8}(undef, ma + mb, length(path))
    i = j = 0
    for (c, op) in enumerate(path)
        if op != _Y
            i += 1
            out[1:ma, c] .= view(a, :, i)
        else
            out[1:ma, c] .= _GAP
        end
        if op != _X
            j += 1
            out[ma+1:end, c] .= view(b, :, j)
        else
            out[ma+1:end, c] .= _GAP
        end
    end
    return out
end

# the best path through the alignment of profiles A and B, as a vector of
# _M / _X / _Y steps
function _gotoh(A::_Profile, B::_Profile, sc::_Scoring)
    LA, LB = length(A.occ), length(B.occ)
    Q = sc.S * B.freq                    # 4 × LB: score of each base against B's columns
    ninf = -Inf32
    Mp = fill(ninf, LB + 1); Xp = fill(ninf, LB + 1); Yp = fill(ninf, LB + 1)
    Mc = similar(Mp); Xc = similar(Xp); Yc = similar(Yp)
    # 2 bits per state: the state each of M, X, Y came from
    tb = zeros(UInt8, LA + 1, LB + 1)
    go, ge = sc.open, sc.extend
    for i in 0:LA
        for j in 0:LB
            if i == 0 && j == 0
                Mc[1], Xc[1], Yc[1] = 0.0f0, ninf, ninf
                continue
            end
            t = 0x00
            # M: columns i and j together
            if i > 0 && j > 0
                s = A.freq[1, i] * Q[1, j] + A.freq[2, i] * Q[2, j] + A.freq[3, i] * Q[3, j] + A.freq[4, i] * Q[4, j]
                v, w = _best3(Mp[j], Xp[j], Yp[j])
                Mc[j+1] = v + s
                t |= w
            else
                Mc[j+1] = ninf
            end
            # X: column i of A against a gap
            if i > 0
                e = ge * A.occ[i]
                o = (j == 0 || j == LB) ? 0.0f0 : go * A.occ[i]
                v, w = _best3(Mp[j+1] - o, Xp[j+1], Yp[j+1] - o)
                Xc[j+1] = v - e
                t |= w << 2
            else
                Xc[j+1] = ninf
            end
            # Y: column j of B against a gap
            if j > 0
                e = ge * B.occ[j]
                o = (i == 0 || i == LA) ? 0.0f0 : go * B.occ[j]
                v, w = _best3(Mc[j] - o, Xc[j] - o, Yc[j])
                Yc[j+1] = v - e
                t |= w << 4
            else
                Yc[j+1] = ninf
            end
            tb[i+1, j+1] = t
        end
        Mp, Mc = Mc, Mp
        Xp, Xc = Xc, Xp
        Yp, Yc = Yc, Yp
    end
    _, state = _best3(Mp[LB+1], Xp[LB+1], Yp[LB+1])
    path = UInt8[]
    sizehint!(path, LA + LB)
    i, j = LA, LB
    while i > 0 || j > 0
        push!(path, state)
        t = tb[i+1, j+1]
        if state == _M
            state = t & 0x03
            i -= 1; j -= 1
        elseif state == _X
            state = (t >> 2) & 0x03
            i -= 1
        else
            state = (t >> 4) & 0x03
            j -= 1
        end
    end
    return reverse!(path)
end

# the largest of the scores of M, X and Y and which it was (M first on ties)
@inline function _best3(m, x, y)
    v, w = m, _M
    x > v && ((v, w) = (x, _X))
    y > v && ((v, w) = (y, _Y))
    return v, w
end

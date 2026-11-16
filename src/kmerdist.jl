# MUSCLE's first-iteration distance for nucleotides, "kmer4_6"
# (fastdistnuc.cpp): MAFFT's 6-tuple distance over the groups A, C, G, T
# and "other" (Katoh et al. 2002).
#
# For each pair, c12 is the sum over the distinct 6-tuples of the second
# sequence of min(count in seq 1, count in seq 2) (counts are bytes and
# wrap at 256, as in MUSCLE); c11 and c22 are the same for each sequence
# with itself. The distance is min(3(c11 − c12)/c11, 3(c22 − c12)/c22).

const _NTUPLE = 6^6

function _tuples(s::Vector{UInt8})
    L = [min(Int(_letter(c)), 4) for c in s]
    n = length(L) - 5
    out = Vector{Int32}(undef, max(n, 0))
    for k in 1:n
        t = 0
        for m in 0:5
            t = 6t + L[k+m]
        end
        out[k] = t
    end
    return out
end

function _count!(cnt::Vector{UInt8}, tup)
    fill!(cnt, 0x00)
    for t in tup
        cnt[t+1] += 0x01
    end
    return cnt
end

function _kmerdist(seqs::Vector{Vector{UInt8}}; threaded::Bool=true, prog=_NOPROGRESS)
    n = length(seqs)
    tups = [_tuples(s) for s in seqs]
    ok = [length(s) >= 5 for s in seqs]
    common = zeros(Int, n, n)
    rowwork = function (i)
        ok[i] || return _next!(prog)
        c1 = _count!(zeros(UInt8, _NTUPLE), tups[i])
        c2 = zeros(UInt8, _NTUPLE)
        for j in 1:i
            ok[j] || continue
            _count!(c2, tups[j])
            s = 0
            for t in tups[j]
                s += min(c1[t+1], c2[t+1])
                c2[t+1] = 0x00
            end
            common[i, j] = common[j, i] = s
        end
        _next!(prog)
    end
    if threaded && Threads.nthreads() > 1 && n > 32
        Threads.@threads :dynamic for i in n:-1:1
            rowwork(i)
        end
    else
        foreach(rowwork, 1:n)
    end
    D = zeros(Float32, n, n)
    for i in 1:n
        c11 = common[i, i] == 0 ? 1.0 : Float64(common[i, i])
        for j in 1:i-1
            c22 = common[j, j] == 0 ? 1.0 : Float64(common[j, j])
            d = min(3.0 * (c11 - common[i, j]) / c11, 3.0 * (c22 - common[i, j]) / c22)
            D[i, j] = D[j, i] = Float32(d)
        end
    end
    return D
end

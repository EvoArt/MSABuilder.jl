# The guide tree: k-mer distances between the unaligned sequences, then
# UPGMA. Only the merge order is needed by the progressive alignment.

# Distance from shared k-mers (as MAFFT's 6-mer distance):
# 1 - Σ min(count_x, count_y) / (k-mers in the shorter sequence).
function _kmerdist(seqs::Vector{Vector{UInt8}}, k::Int; threaded::Bool=true, prog=_NOPROGRESS)
    n = length(seqs)
    K = 4^k
    counts = zeros(UInt8, K, n)
    kmers = Vector{Vector{Int32}}(undef, n)   # distinct k-mers of each sequence
    total = zeros(Int, n)
    for (j, s) in enumerate(seqs)
        c = view(counts, :, j)
        u = Int32[]
        for x in _kmers(s, k)
            c[x+1] == 0 && push!(u, x)
            c[x+1] < typemax(UInt8) && (c[x+1] += 0x01)
            total[j] += 1
        end
        kmers[j] = u
    end
    D = zeros(Float32, n, n)
    rowwork = function (i)
        ci = view(counts, :, i)
        for j in i+1:n
            cj = view(counts, :, j)
            shared = 0
            for x in kmers[i]
                shared += min(ci[x+1], cj[x+1])
            end
            m = min(total[i], total[j])
            D[i, j] = D[j, i] = m == 0 ? 1.0f0 : 1.0f0 - Float32(shared / m)
        end
        _next!(prog)
    end
    if threaded && Threads.nthreads() > 1 && n > 64
        Threads.@threads :dynamic for i in 1:n
            rowwork(i)
        end
    else
        foreach(rowwork, 1:n)
    end
    return D
end

# the k-mers of `s` as integers in 0:4^k-1, skipping any containing a
# letter other than A/C/G/T/U
function _kmers(s::Vector{UInt8}, k::Int)
    out = Int32[]
    mask = Int32(4^k - 1)
    x = Int32(0)
    run = 0
    for b in s
        v = _KCODE[b+1]
        if v == 0xff
            run = 0
            x = Int32(0)
            continue
        end
        x = ((x << 2) | Int32(v)) & mask
        run += 1
        run >= k && push!(out, x)
    end
    return out
end

# UPGMA by the nearest-neighbour chain (O(n²) time). Returns the merges as
# pairs of cluster ids: tips are 1:n, the cluster made by merge t is n + t.
function _upgma(D::Matrix{Float32})
    n = size(D, 1)
    D = copy(D)
    sz = ones(Int, n)
    active = trues(n)
    id = collect(1:n)
    merges = Tuple{Int,Int}[]
    chain = Int[]
    left = n
    while left > 1
        isempty(chain) && push!(chain, findfirst(active))
        a = chain[end]
        prev = length(chain) >= 2 ? chain[end-1] : 0
        # nearest active neighbour of a, keeping prev on ties so the chain ends
        b, bd = prev, prev == 0 ? Inf32 : D[prev, a]
        for c in 1:n
            (active[c] && c != a) || continue
            D[c, a] < bd && ((b, bd) = (c, D[c, a]))
        end
        if b != prev
            push!(chain, b)
            continue
        end
        pop!(chain); pop!(chain)
        push!(merges, (id[prev], id[a]))
        for c in 1:n
            (active[c] && c != a && c != prev) || continue
            D[c, prev] = D[prev, c] = (sz[prev] * D[c, prev] + sz[a] * D[c, a]) / (sz[prev] + sz[a])
        end
        sz[prev] += sz[a]
        active[a] = false
        id[prev] = n + length(merges)
        left -= 1
    end
    return merges
end

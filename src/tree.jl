# Guide trees: MUSCLE's UPGMA2 (upgma2.cpp) with its default "UPGMB"
# linkage, and the ClustalW sequence weights computed from the tree
# (clwwt.cpp; Thompson, Higgins & Gibson 1994).

# A rooted binary tree. Tips are nodes 1:n; internal node n + k is made by
# the k-th join, with children left[k], right[k] and branch lengths
# leftlen[k], rightlen[k]. The root is node 2n − 1.
struct _Tree
    n::Int
    left::Vector{Int}
    right::Vector{Int}
    leftlen::Vector{Float32}
    rightlen::Vector{Float32}
end

const _SUEFF = 0.1f0

# UPGMA2 with LINKAGE_Biased: the distance to a new cluster is
# SUEFF·avg(dL, dR) + (1 − SUEFF)·min(dL, dR). Nearest neighbours are
# tracked as in MUSCLE (including keeping a stale minimum when a
# neighbour is merged), so ties and joins come out the same.
function _upgmb(D::Matrix{Float32})
    n = size(D, 1)
    D = copy(D)
    BIG = 1f29
    mindist = fill(BIG, n)
    nn = zeros(Int, n)
    node = collect(1:n)                  # 0: row deleted
    height = zeros(Float32, 2n - 1)
    t = _Tree(n, zeros(Int, n - 1), zeros(Int, n - 1), zeros(Float32, n - 1), zeros(Float32, n - 1))
    for i in 2:n, j in 1:i-1
        d = D[i, j]
        d < mindist[i] && (mindist[i] = d; nn[i] = j)
        d < mindist[j] && (mindist[j] = d; nn[j] = i)
    end
    for k in 1:n-1
        L = R = 0
        best = BIG
        for j in 1:n
            node[j] == 0 && continue
            if mindist[j] < best
                best = mindist[j]
                L = j
                R = nn[j]
            end
        end
        newmin = BIG
        newnn = 0
        for j in 1:n
            (j == L || j == R || node[j] == 0) && continue
            dL, dR = D[L, j], D[R, j]
            d = _SUEFF * ((dL + dR) / 2) + (1 - _SUEFF) * min(dL, dR)
            nn[j] == R && (nn[j] = L)
            D[L, j] = D[j, L] = d
            d < newmin && (newmin = d; newnn = j)
        end
        h = D[L, R] / 2
        t.left[k], t.right[k] = node[L], node[R]
        t.leftlen[k] = h - height[node[L]]
        t.rightlen[k] = h - height[node[R]]
        height[n+k] = h
        node[L] = n + k
        nn[L] = newnn
        mindist[L] = newmin
        node[R] = 0
    end
    return t
end

# ClustalW weights: a branch's strength is its length over the number of
# tips below it; a tip's weight is the sum of strengths from it to the
# root (1 if that is below 0.0001), normalised to sum to 1.
function _clustalw_weights(t::_Tree)
    n = t.n
    n == 1 && return [1.0f0]
    n == 2 && return [0.5f0, 0.5f0]
    nnode = 2n - 1
    parent = zeros(Int, nnode)
    len = zeros(Float32, nnode)
    leaves = ones(Int, nnode)
    for k in 1:n-1
        v = n + k
        parent[t.left[k]] = v; len[t.left[k]] = t.leftlen[k]
        parent[t.right[k]] = v; len[t.right[k]] = t.rightlen[k]
        leaves[v] = leaves[t.left[k]] + leaves[t.right[k]]
    end
    strength = [parent[v] == 0 ? 0.0 : Float64(len[v]) / leaves[v] for v in 1:nnode]
    w = zeros(Float32, n)
    for i in 1:n
        s = 0.0
        v = i
        while parent[v] != 0
            s += strength[v]
            v = parent[v]
        end
        w[i] = Float32(s < 0.0001 ? 1.0 : s)
    end
    total = 0.0f0
    for x in w
        total += x
    end
    w ./= total
    return w
end

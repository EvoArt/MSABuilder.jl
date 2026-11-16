# Progressive alignment up a guide tree (MUSCLE's ProgressiveAlign,
# progressivealign.cpp, with AlignTwoMSAs, aligntwomsas.cpp).

struct _Params
    matrix::Matrix{Float32}
    gapopen::Float32
    gapextend::Float32
end

# A group of aligned sequences: letters (rows × columns) and the sequence
# ids of its rows.
struct _Group
    rows::Matrix{UInt8}
    ids::Vector{Int}
end

# the ClustalW weights `w` of the group's sequences, normalised to sum 1
# (SetClustalWWeightsMuscle)
function _groupweights(ids, w::Vector{Float32})
    v = w[ids]
    total = 0.0f0
    for x in v
        total += x
    end
    f = 1.0f0 / total
    return v .* f
end

function _aligngroups(a::_Group, b::_Group, w, p::_Params)
    PA = _profile(a.rows, _groupweights(a.ids, w), p)
    PB = _profile(b.rows, _groupweights(b.ids, w), p)
    path = _nwsmall(PA, PB, p.gapextend)
    return _Group(_mergemsa(a.rows, b.rows, path), vcat(a.ids, b.ids))
end

# align `seqs` up `tree`; the result's rows are in the order of `ids`
function _progressive(seqs::Vector{Vector{UInt8}}, tree::_Tree, p::_Params; prog=_NOPROGRESS)
    n = length(seqs)
    w = _clustalw_weights(tree)
    groups = Vector{_Group}(undef, 2n - 1)
    for i in 1:n
        groups[i] = _Group(reshape(copy(seqs[i]), 1, :), [i])
    end
    for k in 1:n-1
        groups[n+k] = _aligngroups(groups[tree.left[k]], groups[tree.right[k]], w, p)
        _next!(prog)
    end
    return groups[2n-1]
end

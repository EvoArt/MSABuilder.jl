# MUSCLE's refinement stages (iteration 2: re-tree from Kimura distances;
# iterations 3+: tree-dependent restricted partitioning) are not ported yet.
function _refine(g, seqs, tree, p, maxiters; kwargs...)
    throw(ArgumentError("maxiters > 1 (MUSCLE's refinement stages) is not available yet; use maxiters = 1"))
end

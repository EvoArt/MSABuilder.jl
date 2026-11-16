# Sequences evolved down a random tree with substitutions, insertions and
# deletions. Every residue carries the id of its true alignment column, so
# an alignment can be scored against the truth. Used by the tests and by
# benchmark/accuracy.jl.

using Random

const _BASES = ('A', 'C', 'G', 'T')

# `seq` evolved along a branch: each site changes with probability `sub`
# (transitions twice as likely as each transversion); Poisson(`indel`) events
# of geometric length (mean `indellen`), insertions and deletions equally
# likely.
function evolve(rng, seq, sub, indel, nextid; indellen=2.0)
    s = map(seq) do (id, b)
        rand(rng) < sub || return (id, b)
        ts = b == 'A' ? 'G' : b == 'G' ? 'A' : b == 'C' ? 'T' : 'C'
        return (id, rand(rng) < 0.5 ? ts : rand(rng, setdiff(_BASES, (b, ts))))
    end
    for _ in 1:_poisson(rng, indel)
        len = 1
        while rand(rng) > 1 / indellen
            len += 1
        end
        length(s) <= len + 2 && continue
        p = rand(rng, 1:length(s)-len+1)
        if rand(rng, Bool)
            deleteat!(s, p:p+len-1)
        else
            for _ in 1:len
                insert!(s, p, (nextid[] += 1, rand(rng, _BASES)))
            end
        end
    end
    return s
end

function _poisson(rng, λ)
    k, p, L = 0, 1.0, exp(-λ)
    while (p *= rand(rng)) > L
        k += 1
    end
    return k
end

"""
    simulate(rng; n=20, L=250, sub=0.1, indel=3.0) -> Vector{Vector{Tuple{Int,Char}}}

`n` sequences from a random root of length `L` evolved down a random binary
tree (random splits, branch lengths uniform in 0.5–1.5, scaled so that a
typical root-to-tip path has length 1). Along that path a site changes with
probability about `sub`, and `indel` insertions or deletions happen.
"""
function simulate(rng; n=20, L=250, sub=0.1, indel=3.0)
    nextid = Ref(L)
    scale = 1 / max(1.0, 2log(n))
    root = [(i, rand(rng, _BASES)) for i in 1:L]
    out = Vector{Vector{Tuple{Int,Char}}}()
    function down(seq, k)
        if k == 1
            push!(out, seq)
            return
        end
        left = rand(rng, 1:k-1)
        for kk in (left, k - left)
            t = (0.5 + rand(rng)) * scale
            down(evolve(rng, seq, sub * t, indel * t, nextid), kk)
        end
    end
    down(root, n)
    return out
end

sequences_of(truth) = [String(last.(s)) for s in truth]

# columns of each residue of each row of an aligned set of strings
function _residuecols(rows)
    map(rows) do r
        cols = Int[]
        for (c, x) in enumerate(r)
            x == '-' || push!(cols, c)
        end
        cols
    end
end

"""
    accuracy(rows, truth) -> (sp, precision, tc)

`sp`: share of the truly homologous residue pairs that the alignment puts
in the same column (sum-of-pairs). `precision`: share of the residue pairs
the alignment puts together that are truly homologous. `tc`: share of the
true columns (with two or more residues) reproduced exactly.
"""
function accuracy(rows::AbstractVector{<:AbstractString}, truth)
    cols = _residuecols(rows)
    ids = [first.(s) for s in truth]
    n = length(truth)
    # true column id -> (row, alignment column) of its residues
    bytrue = Dict{Int,Vector{Tuple{Int,Int}}}()
    byaln = Dict{Int,Vector{Tuple{Int,Int}}}()
    for r in 1:n, (k, id) in enumerate(ids[r])
        push!(get!(bytrue, id, Tuple{Int,Int}[]), (r, cols[r][k]))
        push!(get!(byaln, cols[r][k], Tuple{Int,Int}[]), (r, id))
    end
    tp = truepairs = alnpairs = 0
    tcok = tctot = 0
    for v in values(bytrue)
        m = length(v)
        m < 2 && continue
        truepairs += m * (m - 1) ÷ 2
        counts = Dict{Int,Int}()
        for (_, c) in v
            counts[c] = get(counts, c, 0) + 1
        end
        tp += sum(x -> x * (x - 1) ÷ 2, values(counts))
        tctot += 1
        # exact: all residues in one alignment column, and nothing else there
        if length(counts) == 1 && length(byaln[first(keys(counts))]) == m
            tcok += 1
        end
    end
    for v in values(byaln)
        m = length(v)
        alnpairs += m * (m - 1) ÷ 2
    end
    return (sp=tp / max(1, truepairs), precision=tp / max(1, alnpairs), tc=tcok / max(1, tctot))
end

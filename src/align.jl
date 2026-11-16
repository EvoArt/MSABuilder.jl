# The MSA type and `align`.

"""
    MSA

A multiple sequence alignment, as returned by [`align`](@ref). `msa[i]` is
the `i`th aligned sequence as a `String` (gaps are `-`), `msa["name"]` looks
one up by name, and iterating gives `name => sequence` pairs. `msa.names`
holds the names, `size(msa)` is (sequences, columns). Write it with
[`writemsa`](@ref).

With PhyloBuilder loaded, an `MSA` can be passed straight to its distance
and tree functions (`nj(msa)`, `distances(msa)`, …) or converted with
`Alignment(msa)`.
"""
struct MSA
    names::Vector{String}
    data::Matrix{UInt8}            # sequences × columns, '-' for gaps
    scoring::NamedTuple
end

Base.length(m::MSA) = size(m.data, 1)
Base.size(m::MSA) = size(m.data)
Base.size(m::MSA, d) = size(m.data, d)
Base.getindex(m::MSA, i::Integer) = String(m.data[i, :])
function Base.getindex(m::MSA, s::AbstractString)
    i = findfirst(==(s), m.names)
    i === nothing && throw(KeyError(s))
    return m[i]
end
Base.iterate(m::MSA, i=1) = i > length(m) ? nothing : (m.names[i] => m[i], i + 1)
Base.eltype(::Type{MSA}) = Pair{String,String}
Base.:(==)(a::MSA, b::MSA) = a.names == b.names && a.data == b.data

"""
    align(seqs; names, k=6, match=5, mismatch=-4, gapopen=10, gapextend=1,
          threaded=true, progress) -> MSA
    align(names, seqs; kwargs...)
    align(["name" => seq, …]; kwargs...)
    align(path; kwargs...)                  # a FASTA file

Align nucleotide sequences progressively: a guide tree is built by UPGMA
from `k`-mer distances, then sequences and groups are joined up the tree by
profile alignment with affine gap penalties (a gap of length `n` costs
`gapopen + n * gapextend`; gaps at the ends of the alignment cost
`n * gapextend`).

Sequences can be strings or anything whose `string` is the sequence (e.g.
BioSequences' `LongDNA`). Gaps in the input are removed; IUPAC ambiguity
codes are kept and scored by the bases they allow; `U` is read as `T`.
Names default to `seq1`, `seq2`, …; rows of the result are in input order.

With DenoiseDNA loaded, `align(seqtab)` aligns the ASVs of a
`SequenceTable` (named `ASV1`, `ASV2`, … in column order).
"""
function align(seqs::AbstractVector; names=nothing, kwargs...)
    names = names === nothing ? ["seq$i" for i in eachindex(seqs)] : names
    return align(names, seqs; kwargs...)
end
align(pairs::AbstractVector{<:Pair}; kwargs...) = align(first.(pairs), last.(pairs); kwargs...)
align(path::AbstractString; kwargs...) = align(readsequences(path); kwargs...)

function align(names::AbstractVector, seqs::AbstractVector; k::Integer=6,
               match::Real=5, mismatch::Real=-4, gapopen::Real=10, gapextend::Real=1,
               threaded::Bool=true, progress::Bool=_progress_default())
    length(names) == length(seqs) || throw(ArgumentError("$(length(names)) names for $(length(seqs)) sequences"))
    isempty(seqs) && throw(ArgumentError("no sequences"))
    names = String[string(x) for x in names]
    allunique(names) || throw(ArgumentError("sequence names must be unique"))
    1 <= k <= 12 || throw(ArgumentError("k must be in 1:12"))
    sc = _Scoring(match, mismatch, gapopen, gapextend)
    scoring = (; match, mismatch, gapopen, gapextend, k)
    enc = [_clean(s) for s in seqs]
    n = length(enc)
    n == 1 && return MSA(names, reshape(enc[1], 1, :), scoring)

    bars = progress && n >= 50
    p = _progress(n, "Distances"; unit="rows", enabled=bars)
    D = _kmerdist(enc, Int(k); threaded, prog=p)
    _finish!(p)
    merges = _upgma(D)

    p = _progress(length(merges), "Aligning"; unit="merges", enabled=bars)
    groups = Vector{Matrix{UInt8}}(undef, n + length(merges))
    members = Vector{Vector{Int}}(undef, n + length(merges))
    for i in 1:n
        groups[i] = reshape(enc[i], 1, :)
        members[i] = [i]
    end
    for (t, (a, b)) in enumerate(merges)
        groups[n+t] = _alignprofiles(groups[a], groups[b], sc)
        members[n+t] = vcat(members[a], members[b])
        groups[a] = groups[b] = Matrix{UInt8}(undef, 0, 0)    # free them
        _next!(p)
    end
    _finish!(p)
    data = similar(groups[end])
    data[members[end], :] = groups[end]
    return MSA(names, data, scoring)
end

"""
    writemsa(path, msa; width=60)
    writemsa(io, msa; width=60)

Write `msa` as FASTA, wrapping sequences at `width` characters (`0`: no
wrapping).
"""
function writemsa(io::IO, m::MSA; width::Integer=60)
    for (name, s) in m
        println(io, '>', name)
        if width <= 0
            println(io, s)
        else
            for k in 1:width:length(s)
                println(io, s[k:min(end, k + width - 1)])
            end
        end
    end
end
writemsa(path::AbstractString, m::MSA; kwargs...) = (open(io -> writemsa(io, m; kwargs...), path, "w"); path)

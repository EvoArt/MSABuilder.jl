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
    params::NamedTuple             # MUSCLE options used
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
    align(seqs; names, maxiters=1, gapopen=-400, gapextend=0, threaded=true, progress) -> MSA
    align(names, seqs; kwargs...)
    align(["name" => seq, …]; kwargs...)
    align(path; kwargs...)                  # a FASTA file

Align nucleotide sequences with MUSCLE 3.8.31 (Edgar 2004), ported to
Julia: k-mer distances (MUSCLE's `kmer4_6`), a UPGMB guide tree, ClustalW
sequence weights and progressive profile–profile alignment, with MUSCLE's
DNA scoring (`-gapopen -400`, `-gapextend 0`, half-penalised terminal gaps).
This is MUSCLE's first iteration (`-maxiters 1`); the refinement
iterations are not ported yet, so `maxiters` must be 1.

Sequences can be strings or anything whose `string` is the sequence (e.g.
BioSequences' `LongDNA`). Gaps in the input are removed, `U` is read as
`T`, and characters other than IUPAC nucleotide codes become `N`. Names
default to `seq1`, `seq2`, …; rows of the result are in input order.

With DenoiseDNA loaded, `align(seqtab)` aligns the ASVs of a
`SequenceTable` (named `ASV1`, `ASV2`, … in column order).
"""
function align(seqs::AbstractVector; names=nothing, kwargs...)
    names = names === nothing ? ["seq$i" for i in eachindex(seqs)] : names
    return align(names, seqs; kwargs...)
end
align(pairs::AbstractVector{<:Pair}; kwargs...) = align(first.(pairs), last.(pairs); kwargs...)
align(path::AbstractString; kwargs...) = align(readsequences(path); kwargs...)

function align(names::AbstractVector, seqs::AbstractVector; maxiters::Integer=1,
               gapopen::Real=-400, gapextend::Real=0,
               threaded::Bool=true, progress::Bool=_progress_default())
    length(names) == length(seqs) || throw(ArgumentError("$(length(names)) names for $(length(seqs)) sequences"))
    isempty(seqs) && throw(ArgumentError("no sequences"))
    names = String[string(x) for x in names]
    allunique(names) || throw(ArgumentError("sequence names must be unique"))
    maxiters >= 1 || throw(ArgumentError("maxiters must be ≥ 1"))
    p = _Params(_NUC_SP, gapopen, gapextend)
    params = (; maxiters, gapopen, gapextend)
    enc = [_clean(s) for s in seqs]
    n = length(enc)
    n == 1 && return MSA(names, reshape(enc[1], 1, :), params)

    bars = progress && n >= 50
    pr = _progress(n, "Distances"; unit="rows", enabled=bars)
    tree = _upgmb(_kmerdist(enc; threaded, prog=pr))
    _finish!(pr)
    pr = _progress(n - 1, "Aligning"; unit="joins", enabled=bars)
    g = _progressive(enc, tree, p; prog=pr)
    _finish!(pr)
    if maxiters >= 2 && n > 2
        g = _refine(g, enc, tree, p, maxiters; threaded, bars)
    end
    data = similar(g.rows)
    data[g.ids, :] = g.rows
    return MSA(names, data, params)
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

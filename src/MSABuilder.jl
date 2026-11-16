"""
    MSABuilder

Multiple sequence alignment of nucleotide sequences in pure Julia:
progressive alignment along a k-mer UPGMA guide tree, with profile–profile
Gotoh alignment and affine gap penalties.
"""
module MSABuilder

using Printf: Printf, @sprintf

export MSA, align, readsequences, writemsa

include("style.jl")
include("progress.jl")
include("sequences.jl")
include("guidetree.jl")
include("profile.jl")
include("align.jl")
include("show.jl")

end

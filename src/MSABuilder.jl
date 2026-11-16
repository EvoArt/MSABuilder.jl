"""
    MSABuilder

Multiple sequence alignment of nucleotide sequences: a Julia port of
MUSCLE 3.8.31 (Edgar 2004).
"""
module MSABuilder

using Printf: Printf, @sprintf

export MSA, align, readsequences, writemsa

include("style.jl")
include("progress.jl")
include("alphabet.jl")
include("sequences.jl")
include("kmerdist.jl")
include("tree.jl")
include("profile.jl")
include("progressive.jl")
include("refine.jl")
include("align.jl")
include("show.jl")

end

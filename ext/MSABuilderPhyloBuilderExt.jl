# An MSA can be used wherever PhyloBuilder takes an Alignment.
module MSABuilderPhyloBuilderExt

using MSABuilder: MSA
using PhyloBuilder: PhyloBuilder, Alignment

PhyloBuilder.Alignment(m::MSA) = Alignment(m.names, [m[i] for i in 1:length(m)])

for f in (:distances, :upgma, :wpgma, :nj, :bionj, :basefreqs)
    @eval PhyloBuilder.$f(m::MSA; kwargs...) = PhyloBuilder.$f(Alignment(m); kwargs...)
end
PhyloBuilder.ntaxa(m::MSA) = length(m)
PhyloBuilder.nsites(m::MSA) = size(m, 2)
PhyloBuilder.taxa(m::MSA) = m.names
PhyloBuilder.sequences(m::MSA) = [m[i] for i in 1:length(m)]
PhyloBuilder.writefasta(x, m::MSA; kwargs...) = PhyloBuilder.writefasta(x, Alignment(m); kwargs...)
PhyloBuilder.writephylip(x, m::MSA) = PhyloBuilder.writephylip(x, Alignment(m))

end

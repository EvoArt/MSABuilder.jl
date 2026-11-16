# Align the ASVs of a DenoiseDNA result directly.
module MSABuilderDenoiseDNAExt

using MSABuilder: MSABuilder, align
using DenoiseDNA: SequenceTable, DadaResult

_asvnames(n) = ["ASV$i" for i in 1:n]

MSABuilder.align(st::SequenceTable; names=_asvnames(length(st.sequences)), kwargs...) =
    align(names, st.sequences; kwargs...)
MSABuilder.align(r::DadaResult; names=_asvnames(length(r.sequence)), kwargs...) =
    align(names, r.sequence; kwargs...)

end

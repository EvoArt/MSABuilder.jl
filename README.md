# MSABuilder.jl

[![CI](https://github.com/EvoArt/MSABuilder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/MSABuilder.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://EvoArt.github.io/MSABuilder.jl/dev/)

Multiple sequence alignment of nucleotide sequences: a Julia port of
MUSCLE 3.8.31 (Edgar 2004, *Nucleic Acids Res.* 32:1792). Output matches
MUSCLE 3.8.31 with `-maxiters 1` on the test sets in `test/fixtures`.

Under development: MUSCLE's refinement iterations (`maxiters > 1`) are not
ported yet.

```julia
using Pkg
Pkg.add("MSABuilder")
```

```julia
using MSABuilder

msa = align("asvs.fasta")              # or a vector of sequences, or name => seq pairs
writemsa("asvs_aligned.fasta", msa)
```

With [DenoiseDNA](https://github.com/EvoArt/DenoiseDNA.jl),
[PhyloBuilder](https://github.com/EvoArt/PhyloBuilder.jl) and
[PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl):

```julia
using DenoiseDNA, MSABuilder, PhyloBuilder, PlotPhyloTree

asvs = seqtab_nochim.sequences
tree = align(asvs; names=["ASV$i" for i in eachindex(asvs)]) |> nj
savefig("asv-tree.svg", treeplot(midpoint(tree)))
```

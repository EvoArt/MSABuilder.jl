# MSABuilder.jl

[![CI](https://github.com/EvoArt/MSABuilder.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/EvoArt/MSABuilder.jl/actions/workflows/CI.yml)

Multiple sequence alignment of nucleotide sequences in pure Julia.
Progressive alignment along a k-mer UPGMA guide tree, with profile–profile
alignment and affine gap penalties. No external binaries, so it runs on
Windows, macOS and Linux.

Under development.

```julia
using Pkg
Pkg.add(url="https://github.com/EvoArt/MSABuilder.jl")
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

tree = seqtab_nochim |> align |> nj        # ASVs named ASV1, ASV2, …
savefig("asv-tree.svg", treeplot(midpoint(tree)))
```

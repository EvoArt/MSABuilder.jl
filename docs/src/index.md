# MSABuilder.jl

Multiple sequence alignment of nucleotide sequences in pure Julia.
MSABuilder is a port of MUSCLE 3.8.31 (Edgar 2004): k-mer distances, a
UPGMB guide tree, ClustalW sequence weights and progressive profile–profile
alignment with MUSCLE's DNA scoring. Its output matches MUSCLE 3.8.31 run
with `-maxiters 1`. The refinement iterations (`maxiters > 1`) are not
ported yet.

## Installation

```julia
using Pkg
Pkg.add("MSABuilder")
```

## Usage

```julia
using MSABuilder

msa = align("asvs.fasta")       # a FASTA file, a vector of sequences or name => seq pairs
msa["ASV1"]                     # aligned sequence, gaps as '-'
writemsa("asvs_aligned.fasta", msa)
```

## From reads to a tree

With [DenoiseDNA](https://github.com/EvoArt/DenoiseDNA.jl),
[PhyloBuilder](https://github.com/EvoArt/PhyloBuilder.jl) and
[PlotPhyloTree](https://github.com/EvoArt/PlotPhyloTree.jl), `align` takes a
sequence table and its result goes straight into the tree builders:

```julia
using DenoiseDNA, MSABuilder, PhyloBuilder, PlotPhyloTree

nochim = removeBimeraDenovo(seqtab; method="consensus")
tree   = nochim |> align |> nj      # ASVs are named ASV1, ASV2, … in column order
savefig("asv-tree.svg", treeplot(midpoint(tree)))
```

## Citation

If you use this package, please cite MUSCLE:

Edgar RC (2004). MUSCLE: multiple sequence alignment with high accuracy and
high throughput. *Nucleic Acids Research* 32, 1792–1797.
[doi:10.1093/nar/gkh340](https://doi.org/10.1093/nar/gkh340)

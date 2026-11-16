# MSABuilder.jl

Multiple sequence alignment of nucleotide sequences: a Julia port of
MUSCLE 3.8.31 (Edgar 2004, *Nucleic Acids Res.* 32:1792). No external
binaries, so it runs on Windows, macOS and Linux. Output matches MUSCLE
3.8.31 with `-maxiters 1`; MUSCLE's refinement iterations are not ported
yet.

## Installation

```julia
using Pkg
Pkg.add("MSABuilder")
```

## Example

```@example main
using MSABuilder

msa = align(pkgdir(MSABuilder, "test", "fixtures", "small.fa"))
```

```julia
writemsa("aligned.fasta", msa)
```

An `MSA` is a collection of `name => sequence` pairs, so
[PhyloBuilder](https://github.com/EvoArt/PhyloBuilder.jl) builds trees from
it directly, and with [DenoiseDNA](https://github.com/EvoArt/DenoiseDNA.jl)
loaded `align` takes a sequence table:

```julia
using DenoiseDNA, MSABuilder, PhyloBuilder

tree = seqtab_nochim |> align |> nj        # ASVs named ASV1, ASV2, …
```

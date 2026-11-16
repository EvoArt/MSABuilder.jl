# Writes the unaligned inputs for the MUSCLE parity tests:
#   julia test/fixtures/make_inputs.jl && Rscript test/fixtures/make_fixtures.R
using Random
include(joinpath(@__DIR__, "..", "simulate.jl"))

function writefa(path, seqs)
    open(path, "w") do io
        for (i, s) in enumerate(seqs)
            println(io, ">s", i, '\n', s)
        end
    end
end

# scatter IUPAC codes through some sequences
function ambiguous(rng, seqs)
    map(seqs) do s
        c = collect(s)
        for _ in 1:rand(rng, 0:3)
            c[rand(rng, eachindex(c))] = rand(rng, "NRYKMSW")
        end
        String(c)
    end
end

sets = [
    ("tiny",      Xoshiro(1), (n=3, L=30, sub=0.1, indel=1.0)),
    ("small",     Xoshiro(2), (n=8, L=120, sub=0.1, indel=2.0)),
    ("easy",      Xoshiro(3), (n=20, L=250, sub=0.05, indel=2.0)),
    ("medium",    Xoshiro(4), (n=20, L=250, sub=0.15, indel=4.0)),
    ("hard",      Xoshiro(5), (n=20, L=250, sub=0.3, indel=8.0)),
    ("many",      Xoshiro(6), (n=100, L=250, sub=0.15, indel=4.0)),
    ("long",      Xoshiro(7), (n=10, L=1000, sub=0.15, indel=6.0)),
]
for (name, rng, kw) in sets
    writefa(joinpath(@__DIR__, "$name.fa"), sequences_of(simulate(rng; kw...)))
end
writefa(joinpath(@__DIR__, "ambig.fa"), ambiguous(Xoshiro(8), sequences_of(simulate(Xoshiro(9); n=15, L=200))))
writefa(joinpath(@__DIR__, "short.fa"), ["ACGTACGTAC", "ACG", "ACGTTCGTAC", "ACGTACGAAC", "CGTA"])
writefa(joinpath(@__DIR__, "dups.fa"), vcat(fill("ACGTTGCAACGTAGCTAGCTAGGCT", 3), ["ACGTTGCAACGTAGCTAGGCT", "ACGTTGCATCGTAGCTAGCTAGGCT"]))

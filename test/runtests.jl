using MSABuilder
using Random
using Test

ungap(s) = replace(s, '-' => "")

# ---------------------------------------------------------------- simulation
# Sequences evolved down a two-level tree with substitutions, insertions and
# deletions. Each residue carries the id of its true alignment column, so
# the alignment can be scored against the truth.

function evolve(rng, seq, sub, nindel, nextid)
    s = [(id, rand(rng) < sub ? rand(rng, setdiff("ACGT", b)) : b) for (id, b) in seq]
    for _ in 1:nindel
        len = rand(rng, 1:3)
        p = rand(rng, 2:length(s)-len-1)
        if rand(rng, Bool)
            deleteat!(s, p:p+len-1)
        else
            for k in 1:len
                insert!(s, p, (nextid[] += 1, rand(rng, "ACGT")))
            end
        end
    end
    return s
end

function simulate(rng; L=300, clades=4, leaves=5)
    nextid = Ref(L)
    root = [(i, rand(rng, "ACGT")) for i in 1:L]
    seqs = Vector{Vector{Tuple{Int,Char}}}()
    for _ in 1:clades
        anc = evolve(rng, root, 0.08, 1, nextid)
        for _ in 1:leaves
            push!(seqs, evolve(rng, anc, 0.03, 1, nextid))
        end
    end
    return seqs
end

# share of the truly homologous residue pairs that are in the same column
function sp_score(m::MSA, truth)
    col = map(1:length(m)) do r
        cols = Dict{Int,Int}()
        k = 0
        for (c, x) in enumerate(m[r])
            x == '-' && continue
            k += 1
            cols[truth[r][k][1]] = c
        end
        cols
    end
    hit = tot = 0
    for a in eachindex(truth), b in a+1:length(truth), (id, ca) in col[a]
        cb = get(col[b], id, nothing)
        cb === nothing && continue
        tot += 1
        hit += ca == cb
    end
    return hit / tot
end

@testset "MSABuilder" begin
    @testset "basics" begin
        seqs = ["ACGTACGTTTGACCA", "ACGTACGTGACCA", "ACGAACGTTTGACCA"]
        m = align(seqs; progress=false)
        @test length(m) == 3
        @test m.names == ["seq1", "seq2", "seq3"]
        @test ungap.(collect(last.(m))) == seqs
        @test allequal(length.(last.(m)))
        @test m["seq2"] == m[2]
        @test_throws KeyError m["nope"]

        # identical sequences: no gaps
        m = align(fill("ACGTTGCA", 4); progress=false)
        @test all(m[i] == "ACGTTGCA" for i in 1:4)

        # one sequence; names; pairs; lowercase and gaps in the input
        @test align(["acg-t"]; progress=false)[1] == "ACGT"
        @test align(["a", "b"], ["ACGT", "ACGT"]; progress=false).names == ["a", "b"]
        @test align(["x" => "ACGT", "y" => "ACT"]; progress=false).names == ["x", "y"]
        @test_throws ArgumentError align(String[])
        @test_throws ArgumentError align(["a", "a"], ["ACGT", "ACGT"])
        @test_throws ArgumentError align(["a"], ["ACGT", "ACGT"])
    end

    @testset "gap placement" begin
        # an unambiguous internal deletion
        m = align(["AAACCCGGGTTTACG", "AAACCCTTTACG"]; progress=false)
        @test m[2] == "AAACCC---TTTACG"
        # end gaps pay no opening penalty
        m = align(["ACGTACGTACGGTTCA", "GTACGTACGGTTCA"]; progress=false)
        @test m[2] == "--GTACGTACGGTTCA"
        # ambiguity codes are kept and matched to the bases they allow
        m = align(["ACGTRCGTAC", "ACGTACGTAC"]; progress=false)
        @test m[1] == "ACGTRCGTAC" && m[2] == "ACGTACGTAC"
        @test align(["ACGU", "ACGT"]; progress=false)[1] == "ACGU"
    end

    @testset "accuracy on simulated sequences" begin
        for seed in 1:3
            truth = simulate(Xoshiro(seed))
            m = align([String(last.(s)) for s in truth]; progress=false)
            @test ungap.(collect(last.(m))) == [String(last.(s)) for s in truth]
            @test sp_score(m, truth) > 0.95
        end
    end

    @testset "guide tree" begin
        # UPGMA merges the two closest pairs first
        D = Float32[0 1 4 4; 1 0 4 4; 4 4 0 2; 4 4 2 0]
        mg = MSABuilder._upgma(D)
        @test Set(Set.(mg[1:2])) == Set([Set((1, 2)), Set((3, 4))])
        @test Set(mg[3]) == Set((5, 6))
        # k-mers skip ambiguous letters
        @test length(MSABuilder._kmers(collect(codeunits("ACGTNACGT")), 3)) == 4
    end

    @testset "I/O and display" begin
        m = align(["ACGTACGTTTGACCA", "ACGTACGTGACCA"]; names=["a", "b"], progress=false)
        path = tempname() * ".fasta"
        writemsa(path, m; width=5)
        @test readsequences(path) == ["a" => "ACGTACGTTTGACCA", "b" => "ACGTACGTGACCA"]
        @test align(path; progress=false) == m
        rm(path)
        s = sprint(show, MIME"text/plain"(), m)
        @test occursin("MSA · 2 sequences × 15 columns", s)
        @test occursin(m[2], s)
        @test sprint(show, m) == "MSA(2 sequences × 15 columns)"
    end

    # the extensions, when the packages are installed
    if Base.find_package("PhyloBuilder") !== nothing
        @testset "PhyloBuilder" begin
            @eval using PhyloBuilder
            m = align(["ACGTACGTTTGACCA", "ACGTACGTGACCA", "ACGAACGTTTGACCA", "ACGAACGTTTGACCT"]; progress=false)
            @test Alignment(m) isa Alignment
            @test ntaxa(m) == 4 && nsites(m) == 15
            @test distances(m) == distances(Alignment(m))
            @test nj(m; progress=false) isa PhyloBuilder.PhyloTree
        end
    end
end

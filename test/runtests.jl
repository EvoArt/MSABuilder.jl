using MSABuilder
using Random
using Test

const M = MSABuilder
const FIX = joinpath(@__DIR__, "fixtures")
include("simulate.jl")

ungap(s) = replace(s, '-' => "")
alignedrows(path) = [String(strip(l)) for l in eachline(path) if !startswith(l, '>')]

@testset "MSABuilder" begin
    @testset "parity with MUSCLE 3.8.31 (maxiters = 1)" begin
        # references from test/fixtures/make_fixtures.R (Bioconductor msa)
        for name in ("tiny", "small", "short", "dups", "ambig", "easy", "medium", "hard", "many", "long")
            m = align(joinpath(FIX, "$name.fa"); progress=false)
            @test collect(last.(m)) == alignedrows(joinpath(FIX, "$(name)_it1.fa"))
        end
    end

    @testset "basics" begin
        seqs = ["ACGTACGTTTGACCA", "ACGTACGTGACCA", "ACGAACGTTTGACCA"]
        m = align(seqs; progress=false)
        @test length(m) == 3
        @test m.names == ["seq1", "seq2", "seq3"]
        @test ungap.(collect(last.(m))) == seqs
        @test allequal(length.(last.(m)))
        @test m["seq2"] == m[2]
        @test_throws KeyError m["nope"]
        m4 = align(fill("ACGTTGCA", 4); progress=false)
        @test all(m4[i] == "ACGTTGCA" for i in 1:4)
        @test align(["acg-t"]; progress=false)[1] == "ACGT"
        @test align(["a", "b"], ["ACGT", "ACGT"]; progress=false).names == ["a", "b"]
        @test align(["x" => "ACGT", "y" => "ACT"]; progress=false).names == ["x", "y"]
        @test align(["ACGTRCGTAC", "ACGTACGTAC"]; progress=false)[1] == "ACGTRCGTAC"
        @test align(["ACGU", "ACGT"]; progress=false)[1] == "ACGT"           # U read as T
        @test align(["AC!GT", "ACNGT"]; progress=false)[1] == "ACNGT"        # unknown → N
        @test_throws ArgumentError align(String[])
        @test_throws ArgumentError align(["a", "a"], ["ACGT", "ACGT"])
        @test_throws ArgumentError align(["a"], ["ACGT", "ACGT"])
        @test_throws ArgumentError align(["ACGT", "ACGA", "ACGG"]; maxiters=2, progress=false)
        @test_throws ArgumentError align(["ACGT", "ACGA"]; maxiters=0)
    end

    @testset "k-mer distance" begin
        s = M._clean.(["ACGTACGTAC", "ACGTACGTAC", "TTTTTTTTTT", "ACG"])
        D = M._kmerdist(s; threaded=false)
        @test D[1, 2] == 0
        @test D[1, 3] == 3                     # no shared 6-tuples
        @test D == D'
        @test length(M._tuples(M._clean("ACGTAC"))) == 1
        @test length(M._tuples(M._clean("ACGTNACGT"))) == 4    # N is in the "other" group
        @test M._kmerdist(s; threaded=true) == D
    end

    @testset "UPGMB guide tree and ClustalW weights" begin
        D = Float32[0 1 4 4; 1 0 4 4; 4 4 0 2; 4 4 2 0]
        t = M._upgmb(D)
        @test (t.left[1], t.right[1]) == (1, 2)
        @test (t.left[2], t.right[2]) == (3, 4)
        @test Set((t.left[3], t.right[3])) == Set((5, 6))
        @test t.leftlen[1] == 0.5f0 && t.rightlen[1] == 0.5f0
        w = M._clustalw_weights(t)
        @test sum(w) ≈ 1
        @test w[1] == w[2] && w[3] == w[4] && w[1] < w[3]
        @test M._clustalw_weights(M._upgmb(Float32[0 1; 1 0])) == [0.5f0, 0.5f0]
    end

    @testset "profiles and NWSmall" begin
        p = M._Params(M._NUC_SP, -400.0f0, 0.0f0)
        P = M._profile(reshape(collect(codeunits("AC-G")), 1, :), [1.0f0], p)
        @test P.counts[:, 1] == [1, 0, 0, 0] && P.counts[:, 3] == zeros(4)
        @test P.gapopen[3] == 0 && P.gapopen[1] == -200
        # wildcards: X is split between G and A (as in MUSCLE), others add w/20 to each base
        Q = M._profile(reshape(collect(codeunits("XN")), 1, :), [1.0f0], p)
        @test Q.counts[:, 1] == Float32[0.5, 0, 0.5, 0]
        @test Q.counts[:, 2] == fill(0.05f0, 4)
        @test align(["AAACCCGGGTTTACG", "AAACCCTTTACG"]; progress=false)[2] == "AAACCC---TTTACG"
        @test align(["ACGTACGTACGGTTCA", "GTACGTACGGTTCA"]; progress=false)[2] == "--GTACGTACGGTTCA"
    end

    @testset "accuracy on simulated sequences" begin
        for (seed, kw) in ((1, (n=20, sub=0.05, indel=2.0)), (2, (n=20, sub=0.15, indel=4.0)), (3, (n=60, sub=0.15, indel=4.0)))
            truth = simulate(Xoshiro(seed); kw...)
            m = align(sequences_of(truth); progress=false)
            @test ungap.(collect(last.(m))) == sequences_of(truth)
            @test accuracy(collect(last.(m)), truth).sp > 0.9
        end
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
        @test occursin("MUSCLE 3.8.31", s)
        @test occursin(m[2], s)
        @test sprint(show, m) == "MSA(2 sequences × 15 columns)"
    end

    # PhyloBuilder takes an MSA as a collection of name => sequence pairs
    if Base.find_package("PhyloBuilder") !== nothing
        @testset "PhyloBuilder" begin
            @eval using PhyloBuilder
            m = align(["ACGTACGTTTGACCA", "ACGTACGTGACCA", "ACGAACGTTTGACCA", "ACGAACGTTTGACCT"]; progress=false)
            @test Alignment(m) isa Alignment
            @test ntaxa(Alignment(m)) == 4 && nsites(Alignment(m)) == size(m, 2)
            @test distances(m) == distances(Alignment(m))
            @test nj(m; progress=false) isa PhyloBuilder.PhyloTree
        end
    end
end

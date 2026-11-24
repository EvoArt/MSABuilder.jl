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

    @testset "coloured display, long and many sequences" begin
        rng = MersenneTwister(3)
        seqs = [join(rand(rng, "ACGT", 150)) for _ in 1:12]
        m = align(seqs; names=["sequence_with_a_very_long_name_$i" for i in 1:12], progress=false)
        io = IOContext(IOBuffer(), :color => true, :displaysize => (24, 80))
        show(io, MIME"text/plain"(), m)
        s = String(take!(io.io))
        @test occursin("\e[", s)                          # coloured
        @test occursin("…", s)                            # long rows and names shortened
        @test occursin("4 more sequences", s)             # 8 of 12 shown
        io = IOContext(IOBuffer(), :rows => 12)
        show(io, MIME"text/plain"(), m)
        @test !occursin("more sequences", String(take!(io.io)))
        @test M._pct(NaN) == "–" && M._pct(1.0) == "100%" && M._pct(1e-4) == "<0.1%"
        @test M._g(NaN) == "NaN" && M._g(Inf) == "Inf" && M._g(-Inf) == "-Inf" && M._g(0.0) == "0"
        @test M._g(12345.0) == "1.23e+04" && M._g(0.5) == "0.5"
        @test M._fmt(-1234567) == "-1,234,567" && M._fmt(3.0) == "3" && M._fmt(0.25) == "0.25"
        @test M._bin([1.0, NaN, 3.0, 5.0], 2) == [1.0, 4.0]
        @test M._bin([NaN, NaN, 1.0], 2) == [NaN, 1.0] || isequal(M._bin([NaN, NaN, 1.0], 2), [NaN, 1.0])
        @test occursin(' ', M._spark(IOContext(IOBuffer(), :color => true), [0.0, NaN, 1.0], 0, 1))
        @test M._spark(IOBuffer(), [2.0, 2.0], 2, 2) == "▅▅"
    end

    @testset "progress bars" begin
        M.setprogress!(true); @test M._progress_default()
        M.setprogress!(false); @test !M._progress_default()
        M._PROGRESS[] = nothing
        withenv("CI" => "true") do
            @test !M._progress_default()
        end
        withenv("CI" => nothing) do
            @test M._progress_default() isa Bool
        end
        for n in (5, 0), color in (false, true)
            buf = IOBuffer()
            io = IOContext(buf, :color => color, :displaysize => (24, 100))
            p = M._progress(n, "Testing"; unit="items", enabled=true, io=io)
            p.tty = true                                   # draw as on a terminal
            M._render(p)
            p.tlast = 0.0; M._next!(p; info="first")
            p.tlast = 0.0; M._next!(p, 2)
            M._finish!(p; info="done", ok=color)
            M._finish!(p)                                  # a second finish does nothing
            out = String(take!(buf))
            @test occursin("Testing", out) && occursin("done", out) && occursin("\r", out)
        end
        # not a terminal: only the final line, and a partial count on failure
        buf = IOBuffer()
        p = M._progress(3, "Quiet"; enabled=true, io=buf)
        M._next!(p, 2); M._finish!(p; ok=false)
        out = String(take!(buf))
        @test count(==('\n'), out) == 1 && occursin("2", out)
        q = M._progress(3, "Off"; enabled=false, io=devnull)
        @test M._next!(q) === q && M._finish!(q) === q
        @test M._secs(0.5) == "500 ms" && M._secs(5.0) == "5.0 s"
        @test M._secs(125.0) == "2 min 5 s" && M._secs(7300.0) == "2 h 2 min"
        # the bars of a real alignment (50+ sequences), written to a non-terminal stderr
        rng = MersenneTwister(4)
        err = mktemp() do path, f
            redirect_stderr(f) do
                align([join(rand(rng, "ACGT", 30)) for _ in 1:50]; progress=true)
            end
            close(f); read(path, String)
        end
        @test occursin("Distances", err) && occursin("Aligning", err)
    end

    @testset "DenoiseDNA extension" begin
        @eval using DenoiseDNA, DataFrames
        seqs = ["ACGTACGTTTGACCA", "ACGTACGTGACCA"]
        st = DenoiseDNA.SequenceTable(["s1"], seqs, [3 2])
        @test align(st; progress=false).names == ["ASV1", "ASV2"]
        r = DenoiseDNA.DadaResult(seqs, [3, 2], DataFrame(), [1, 2], [0.0, 0.0], zeros(Int, 16, 2),
                                  zeros(0, 0), DataFrame(), nothing, nothing, DenoiseDNA.DadaOptions())
        @test align(r; names=["x", "y"], progress=false).names == ["x", "y"]
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

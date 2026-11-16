# Accuracy and speed on simulated alignments with a known truth.
#
#   julia --project=benchmark benchmark/accuracy.jl [quick|full]
#
# benchmark/Project.toml has MSABuilder (developed from ..) and MUSCLE_jll,
# which is used as a baseline when it loads.

using MSABuilder, Printf, Random
include(joinpath(@__DIR__, "..", "test", "simulate.jl"))

const HAVE_MUSCLE = try
    @eval using MUSCLE_jll
    true
catch
    false
end

function run_muscle(seqs)
    dir = mktempdir()
    inp, out = joinpath(dir, "in.fa"), joinpath(dir, "out.fa")
    open(inp, "w") do io
        for (i, s) in enumerate(seqs)
            println(io, ">s", i, '\n', s)
        end
    end
    t = @elapsed run(pipeline(`$(MUSCLE_jll.muscle()) -align $inp -output $out`; stdout=devnull, stderr=devnull))
    rows = Dict(readsequences(inp) .|> first .=> "")
    name = ""
    for line in eachline(out)
        startswith(line, '>') ? (name = String(strip(line[2:end]))) : (rows[name] *= strip(line))
    end
    return [rows["s$i"] for i in eachindex(seqs)], t
end

function scenario(; n, L, sub, indel, reps, muscle)
    res = Dict{String,Vector{NTuple{4,Float64}}}()
    for rep in 1:reps
        truth = simulate(Xoshiro(rep); n, L, sub, indel)
        seqs = sequences_of(truth)
        t = @elapsed m = align(seqs; progress=false)
        a = accuracy(collect(last.(m)), truth)
        push!(get!(res, "MSABuilder", []), (a.sp, a.precision, a.tc, t))
        if muscle && HAVE_MUSCLE
            rows, t = run_muscle(seqs)
            a = accuracy(rows, truth)
            push!(get!(res, "MUSCLE 5", []), (a.sp, a.precision, a.tc, t))
        end
    end
    for (tool, v) in sort(collect(res))
        mean(i) = sum(x[i] for x in v) / length(v)
        @printf("%5d %5d %5.2f %5.1f  %-11s %6.3f %6.3f %6.3f %8.3f\n", n, L, sub, indel, tool, mean(1), mean(2), mean(3), mean(4))
    end
end

mode = isempty(ARGS) ? "quick" : ARGS[1]
align(sequences_of(simulate(Xoshiro(0); n=10)); progress=false)     # compile
println("    n     L   sub indel  tool            SP   prec     TC   time(s)")
grid = mode == "quick" ?
    [(n=20, L=250, sub=s, indel=i, reps=5, muscle=true) for (s, i) in ((0.05, 2.0), (0.15, 4.0), (0.3, 8.0))] :
    vcat([(n=n, L=250, sub=s, indel=i, reps=3, muscle=n <= 200)
          for n in (20, 100, 200, 500) for (s, i) in ((0.05, 2.0), (0.15, 4.0), (0.3, 8.0))],
         [(n=50, L=1000, sub=0.15, indel=8.0, reps=3, muscle=true)])
for g in grid
    scenario(; g...)
end

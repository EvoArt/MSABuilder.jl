# Reading unaligned sequences.

"""
    readsequences(path) -> Vector{Pair{String,String}}

Read the sequences of a FASTA file as `name => sequence` pairs, with gaps
and whitespace removed. The name is the first word of the header when that
keeps names unique, the whole header otherwise.
"""
function readsequences(path::AbstractString)
    names = String[]
    seqs = String[]
    buf = IOBuffer()
    for line in eachline(path)
        line = rstrip(line)
        if startswith(line, '>')
            isempty(names) || push!(seqs, String(take!(buf)))
            push!(names, String(strip(line[2:end])))
        elseif !isempty(names)
            print(buf, line)
        end
    end
    isempty(names) || push!(seqs, String(take!(buf)))
    isempty(names) && throw(ArgumentError("$path: no FASTA records"))
    ids = [isempty(n) ? n : String(first(split(n))) for n in names]
    names = allunique(ids) ? ids : names
    return [n => String(_clean(s)) for (n, s) in zip(names, seqs)]
end

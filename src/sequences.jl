# Input sequences: reading, cleaning and encoding.

const _GAPCHARS = (UInt8('-'), UInt8('.'))
const _GAP = UInt8('-')

# IUPAC nucleotide codes as fractions of A, C, G, T for scoring. N, gaps
# and unknown letters score 0 against everything.
const _FRAC = let t = zeros(Float32, 4, 256)
    codes = Dict('A' => "A", 'C' => "C", 'G' => "G", 'T' => "T", 'U' => "T",
                 'R' => "AG", 'Y' => "CT", 'S' => "CG", 'W' => "AT", 'K' => "GT", 'M' => "AC",
                 'B' => "CGT", 'D' => "AGT", 'H' => "ACT", 'V' => "ACG")
    for (c, bases) in codes, b in bases
        t[findfirst(b, "ACGT"), UInt8(c)+1] = 1 / length(bases)
    end
    t
end

# 2-bit code of A/C/G/T(U) for k-mers; 0xff otherwise
const _KCODE = let t = fill(0xff, 256)
    for (c, v) in zip("ACGTU", (0, 1, 2, 3, 3))
        t[UInt8(c)+1] = v
    end
    t
end

# uppercase ASCII bytes of `s` without gaps or whitespace
function _clean(s)
    str = string(s)
    out = UInt8[]
    sizehint!(out, ncodeunits(str))
    for b in codeunits(str)
        (b in _GAPCHARS || isspace(Char(b))) && continue
        b < 0x80 || throw(ArgumentError("non-ASCII character in a sequence"))
        push!(out, UInt8(uppercase(Char(b))))
    end
    return out
end

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

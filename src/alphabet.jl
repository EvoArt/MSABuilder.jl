# Nucleotide letters as in MUSCLE 3.8.31 (alpha.h, alpha.cpp): A C G T are
# letters 0–3; the IUPAC ambiguity codes are wildcards. U is read as T.
# Other characters become N, as MUSCLE's FixAlpha does.

const _GAP = UInt8('-')
const _GAPCHARS = (UInt8('-'), UInt8('.'))

const _NX_A, _NX_C, _NX_G, _NX_T = 0, 1, 2, 3
const _NX_X = 14                       # 'X'; see _addcounts! in profile.jl

# letter of each (uppercase) character: 0–3 for A C G T, 4–15 for the
# wildcards M R W S Y K V H D B X N (MUSCLE's NX enum)
const _NXCODE = let t = fill(Int8(15), 256)
    for (i, c) in enumerate("ACGTMRWSYKVHDBXN")
        t[UInt8(c)+1] = i - 1
    end
    t[UInt8('U')+1] = _NX_T
    t
end
@inline _letter(c::UInt8) = _NXCODE[c+1]

# uppercase ASCII letters of `s`, without gaps or whitespace; U → T and
# characters that are not IUPAC nucleotide codes → N
function _clean(s)
    str = string(s)
    out = UInt8[]
    sizehint!(out, ncodeunits(str))
    for b in codeunits(str)
        (b in _GAPCHARS || isspace(Char(b))) && continue
        b < 0x80 || throw(ArgumentError("non-ASCII character in a sequence"))
        c = UInt8(uppercase(Char(b)))
        c = c == UInt8('U') ? UInt8('T') : c
        push!(out, occursin(Char(c), "ACGTMRWSYKVHDBXN") ? c : UInt8('N'))
    end
    return out
end

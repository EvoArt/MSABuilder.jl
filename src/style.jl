# Terminal styling shared by the object displays (show.jl) and the progress
# bars (progress.jl), in the same look as DenoiseDNA.jl and PhyloBuilder.jl.
#
# Colours are xterm-256 indices: they look the same in every modern terminal,
# whereas the 16 named colours are re-themed per terminal. Escape codes are
# written only when the IO has :color, so sprint(), logs and files get plain
# text.

const _PAL = (head=117, key=245, rule=240, dim=242, good=114, warn=221, amber=214, bad=203,
              spark=117, tail=239, hi=255,
              A=114, C=75, G=214, T=203, N=245, gap=239)

_hascolor(io) = get(io, :color, false)::Bool

# `s` in an xterm-256 colour (and/or bold), or plain `s` when the IO has no colour
function _sty(io, s; color=nothing, bold::Bool=false)
    _hascolor(io) || return string(s)
    pre = (bold ? "\e[1m" : "") * (color === nothing ? "" : "\e[38;5;$(color)m")
    return isempty(pre) ? string(s) : string(pre, s, "\e[0m")
end

_termwidth(io) = clamp(displaysize(io)[2], 40, 200)

# 1234567 -> "1,234,567"
function _fmt(n::Integer)
    s = string(abs(n))
    k = length(s)
    b = IOBuffer()
    n < 0 && print(b, '-')
    for (i, c) in enumerate(s)
        print(b, c)
        r = k - i
        r > 0 && r % 3 == 0 && print(b, ',')
    end
    return String(take!(b))
end
_fmt(x::Real) = isinteger(x) && abs(x) < 1e15 ? _fmt(Int(x)) : _g(x)

_pct(x) = isnan(x) ? "–" : x >= 0.9995 ? "100%" : x < 0.001 && x > 0 ? "<0.1%" : @sprintf("%.1f%%", 100x)
_g(x) = isnan(x) ? "NaN" : isinf(x) ? (x > 0 ? "Inf" : "-Inf") :
        x == 0 ? "0" : abs(x) >= 1e4 || abs(x) < 1e-2 ? @sprintf("%.2e", x) : @sprintf("%.3g", x)

# "── Title ───────" to `width`
function _section(io, title, width)
    t = " " * title * " "
    println(io)
    println(io, _sty(io, "──"; color=_PAL.rule), _sty(io, t; color=_PAL.head, bold=true),
            _sty(io, repeat('─', max(2, width - 2 - textwidth(t))); color=_PAL.rule))
end

# the first line of a display: the type name, then what it holds
_title(io, name, rest) = println(io, _sty(io, name; color=_PAL.head, bold=true), _sty(io, " · "; color=_PAL.rule), rest)

const _KEY_W = 12
_key(io, k) = _sty(io, "  " * rpad(k, _KEY_W); color=_PAL.key)
_kv(io, k, v) = println(io, _key(io, k), v)
_note(io, s) = _sty(io, s; color=_PAL.key)

# a nucleotide sequence in colour, shortened to `w` characters around a
# middle ellipsis
function _seq(io, s::AbstractString, w::Int)
    if length(s) > w
        h = (w - 1) ÷ 2
        t = w - 1 - h
        return _bases(io, first(s, h)) * _sty(io, "…"; color=_PAL.rule) * _bases(io, last(s, t))
    end
    return _bases(io, s)
end

function _bases(io, s::AbstractString)
    _hascolor(io) || return String(s)
    b = IOBuffer()
    prev = nothing
    for c in s
        col = c == 'A' ? _PAL.A : c == 'C' ? _PAL.C : c == 'G' ? _PAL.G : c == 'T' || c == 'U' ? _PAL.T :
              c == '-' ? _PAL.gap : _PAL.N
        col == prev || print(b, "\e[38;5;$(col)m")
        print(b, c)
        prev = col
    end
    print(b, "\e[0m")
    return String(take!(b))
end

const _SPARK = ('▁', '▂', '▃', '▄', '▅', '▆', '▇', '█')

# one glyph per value of `v` on the scale `lo..hi`; `color(k, x)` picks each
# glyph's colour (runs of one colour are written as one escape sequence)
function _spark(io, v::AbstractVector{<:Real}, lo::Real, hi::Real; color=(k, x) -> _PAL.spark)
    b = IOBuffer()
    prev = -1
    for (k, x) in enumerate(v)
        if isnan(x)
            g, c = ' ', 0
        else
            f = hi > lo ? clamp((x - lo) / (hi - lo), 0, 1) : 0.5
            g, c = _SPARK[1+round(Int, 7f)], color(k, x)
        end
        if _hascolor(io) && c != prev
            print(b, c == 0 ? "\e[0m" : "\e[38;5;$(c)m")
            prev = c
        end
        print(b, g)
    end
    _hascolor(io) && prev != -1 && print(b, "\e[0m")
    return String(take!(b))
end

# average `v` into at most `w` bins (NaNs ignored)
function _bin(v::AbstractVector{<:Real}, w::Int)
    n = length(v)
    n <= w && return collect(Float64, v)
    out = fill(NaN, w)
    for k in 1:w
        lo = div((k - 1) * n, w) + 1
        hi = div(k * n, w)
        s = 0.0; c = 0
        for i in lo:hi
            isnan(v[i]) && continue
            s += v[i]; c += 1
        end
        c > 0 && (out[k] = s / c)
    end
    return out
end


# rows to show of `n`; `:rows => k` in the IOContext asks for more
_nshow(io, n, k) = min(n, get(io, :rows, k)::Int)
_more(io, n, shown, what) = n > shown && println(io, _note(io, "  … $(_fmt(n - shown)) more $what"))

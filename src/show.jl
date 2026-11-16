# Displays: `show(io, x)` is one line; `show(io, MIME"text/plain"(), x)`
# (what the REPL prints) is a headline, key-value lines and the first rows,
# sized to the terminal. Helpers are in style.jl.

Base.show(io::IO, m::MSA) = print(io, "MSA(", _fmt(length(m)), " sequences × ", _fmt(size(m, 2)), " columns)")

# share of each column's sequences carrying its most common letter
function _conservation(m::MSA)
    n, L = size(m)
    out = zeros(L)
    counts = zeros(Int, 256)
    for c in 1:L
        fill!(counts, 0)
        for r in 1:n
            x = m.data[r, c]
            x == _GAP || (counts[x+1] += 1)
        end
        out[c] = maximum(counts) / n
    end
    return out
end

function Base.show(io::IO, ::MIME"text/plain", m::MSA)
    n, L = size(m)
    _title(io, "MSA", "$(_fmt(n)) sequences × $(_fmt(L)) columns")
    n == 0 && return
    lens = [count(!=(_GAP), view(m.data, r, :)) for r in 1:n]
    lo, hi = extrema(lens)
    _kv(io, "lengths", lo == hi ? "$lo bp" : "$(lo)–$(hi) bp")
    ngap = count(==(_GAP), m.data)
    gapcols = count(c -> any(==(_GAP), view(m.data, :, c)), 1:L)
    _kv(io, "gaps", "$(_pct(ngap / max(1, n * L))) of cells" * _note(io, " · $(_fmt(gapcols)) columns with gaps"))
    cons = _conservation(m)
    _kv(io, "identical", "$(_fmt(count(==(1), cons))) columns" * _note(io, " · $(_pct(count(==(1), cons) / max(1, L)))"))
    W = _termwidth(io)
    if L > 0
        col(k, x) = x >= 0.9 ? _PAL.good : x >= 0.5 ? _PAL.amber : _PAL.bad
        _kv(io, "conserved", _spark(io, _bin(cons, max(10, min(W, 100) - _KEY_W - 4)), 0, 1; color=col))
    end
    s = m.params
    _kv(io, "method", "MUSCLE 3.8.31" * _note(io, " · maxiters $(s.maxiters), gap open $(s.gapopen), extend $(s.gapextend)"))
    k = _nshow(io, n, 8)
    _section(io, "Sequences", min(W, 100))
    nw = min(24, maximum(length, view(m.names, 1:k)))
    for i in 1:k
        nm = length(m.names[i]) > nw ? first(m.names[i], nw - 1) * "…" : m.names[i]
        println(io, "  ", _note(io, rpad(nm, nw)), "  ", _seq(io, m[i], max(12, W - nw - 6)))
    end
    _more(io, n, k, "sequences")
end

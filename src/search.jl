module search

UR = UnitRange{Int}

function search1d(f::Function, x::UR)::Union{UR, Nothing}
    fl = f(x.start)
    x.start == x.stop && return fl ? x : nothing
    fh = f(x.stop)
    if fl
        fh && return x
        return x.start:search1dtf(f, x.start, x.stop)
    elseif fh
        return search1dft(f, x.start, x.stop):x.stop
    end
    return search1dff(f, x.start, x.stop)
end

function search1dtf(f::Function, lo::Int, hi::Int)::Int
    mid = div(lo + hi, 2)
    mid == lo && return lo
    f(mid) ? search1dtf(f, mid, hi) : search1dtf(f, lo, mid)
end

function search1dft(f::Function, lo::Int, hi::Int)::Int
    mid = div(lo + hi, 2)
    mid == lo && return hi
    f(mid) ? search1dft(f, lo, mid) : search1dft(f, mid, hi)
end

function search1dff(f::Function, lo::Int, hi::Int)::Union{UR, Nothing}
    mid = div(lo, hi, 2)
    mid == lo && return nothing
    f(mid) && search1dft(f, lo, mid):search1dtf(f, mid, hi)
    search1dff(lo, mid) || search1dff(mid + 1, hi)
end

function search(f::Function, x::UR, y::UR)::Vector{Tuple{UR,UR}}
    length(x) > length(y) && return reverse.(search(f, y, x))
    y_xlo = search1d(y -> f(x.start, y), y)
    if x.start == x.stop
        y_xlo == nothing && return []
        return [(x, y_xlo)]
    end
    y_xhi = search1d(y -> f(x.stop, y), y)
    if y_xlo == nothing
        y_xhi == nothing && return searchff(f, x, y)
        return searchft(f, x, y, y_xhi)
    elseif y_xhi == nothing
        return searchtf(f, x, y, y_xlo)
    elseif y_xlo.start == y_xhi.start
        y_xlo.stop == y_xhi.stop && return [(x,
        if y_xlo.stop > y_xhi.stop
            corner = search1dtf(x -> f(x, y_xlo.stop), y)
            return [(





    xhi = search1d(
    pred(x.start, y.start)

end

end # module

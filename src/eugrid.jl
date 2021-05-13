module eugrid

export shortest_paths, Grid

import LinearAlgebra

function shortest_paths(diags::BitMatrix)::Tuple{Matrix{Int},Matrix{Set{CartesianIndex{2}}}}
    n, m = size(diags)
    d = (0:n) .+ (0:m)'
    b = similar(d, Set{CartesianIndex{2}})
    b[:, 1] .= (eltype(b)(),)
    b[1, :] .= (b[1, 1],)

    for i = 1:n
        for j = 1:m
            shortest = min(d[i+1, j], d[i, j+1])
            if diags[i, j] && d[i, j] <= shortest
                if d[i, j] == shortest
                    b[i+1, j+1] = intersect(b[i+1, j], b[i, j], b[i, j+1])
                elseif d[i, j] <= shortest
                    shortest = d[i, j]
                    b[i+1, j+1] = union(b[i, j], (CartesianIndex(i, j),))
                end
            elseif d[i+1, j] == d[i, j+1]
                b[i+1, j+1] = intersect(b[i+1, j], b[i, j+1])
            elseif d[i+1, j] == shortest
                b[i+1, j+1] = b[i+1, j]
            else
                b[i+1, j+1] = b[i, j+1]
            end
            d[i+1, j+1] = shortest + 1
        end
    end
    d, b
end

lindices(n::Int) = collect(zip([1:n; fill(n, n - 1)], [fill(0, n); 1:n-1]))
tindices(n::Int) = reverse.(lindices(n))

taxicab(v, w)::Int16 = LinearAlgebra.norm(v .- w, 1)
euclid(v, w)::Float64 = LinearAlgebra.norm(v .- w, 2)

struct Grid
    n::Int
    diags::BitMatrix

    dil::Array{Int16,3}
    dit::Array{Int16,3}

    dg::Matrix{Int}
    de::Matrix{Float64}
    dd::Matrix{Float64}

    function Grid(n::Int)
        iindices = CartesianIndices((n, 0:n-1)) .|> x -> x.I

        dil = [taxicab(inner, outer) for inner in iindices, outer in lindices(n)]
        dit = permutedims(dil, (2, 1, 3))

        dg = [taxicab(l, t) for l in lindices(n), t in tindices(n)]
        de = [euclid(l, t) for l in lindices(n), t in tindices(n)]
        dd = cumsum(cumsum(2 .* (dg .- de) .- 1, dims=1), dims=2)

        new(n, falses(n, n), dil, dit, dg, de, dd)
    end
end

struct Impact
    l::UnitRange{Int}
    t::UnitRange{Int}
end

function score(g::Grid, i::Impact)::Float64
    s = g.dd[i.l.stop, i.t.stop]
    if i.l.start > 1
        s -= g.dd[i.l.start - 1, i.t.stop]
    end
    if i.t.start > 1
        s -= g.dd[i.l.stop, i.t.start - 1]
    end
    if i.l.start > 1 && i.t.start > 1
        s += g.dd[i.l.start - 1, i.t.start - 1]
    end
    s
end

struct Flip
    x::Int
    y::Int
    impacts::Vector{Impact}

    Flip(n::Int, x::Int, y::Int) = new(x, y, [Impact(x:n+y-1,y:n+x-1)])
end

score(g::Grid, f::Flip)::Float64 = sum(score(g, i) for i in f.impacts)

struct Flipper
    g::Grid
    flips::Matrix{Flip}

    Flipper(n::Int) = new(Grid(n), [Flip(n, i.I...) for i in CartesianIndices((n, n))])
end

score(f::Flipper)::Matrix{Float64} = f.flips .|> x->score(f.g, x)

function flip(f::Flipper, x::Int, y::Int)::Nothing

end


    #=
function add(g::Grid, x::Int, y::Int)::Nothing
    @assert !g.diags[x, y]
    g.diags[x, y] = true

    lix = x:g.n_inner+y-1
    dxyl = g.dil[x, y, lix]

    if x > 1 && y < g.n_inner
        sps, bb = shortest_paths(g.diags[x-1:-1:2, y+1:g.n_inner-1])
        dil = view(g.dil, x-1:-1:1, y+1:g.n_inner, lix)
        dil .= min.(dil, reshape(dxyl, 1, 1, :) .+ sps .+ 1)
    end

    tix = y:g.n_inner+x-1
    dxyt = g.dit[x, y, tix]

    if x < g.n_inner && y > 1
        sps, bb = shortest_paths(g.diags[x+1:g.n_inner-1, y-1:-1:2])
        dit = view(g.dit, x+1:g.n_inner, y-1:-1:1, tix)
        dit .= min.(dit, reshape(dxyt, 1, 1, :) .+ sps)
    end

    dg = view(g.dg, lix, tix)
    dg .= min.(dg, reshape(dxyl, :, 1) .+ reshape(dxyt, 1, :) .+ 1)

    return
end


function apply(g::Grid, ::Add, given::Add)::Float64

        Threads.@threads for ij in  CartesianIndices(scores)
            i, j = ij.I
            #delta = view(deltas, i:n_inner+j-1, j:n_inner+i-1)
            #scores[ij] = sum(delta)
            scores[ij] = sum(
                deltas[d] for d in CartesianIndices((i:n_inner+j-1, j:n_inner+i-1)))
        end

=#

end # module

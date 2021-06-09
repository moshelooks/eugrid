module eugrid

import LinearAlgebra

using CircularArrays: CircularArray

#export shortest_paths, Grid

include("antidiagonal.jl")
include("torus.jl")
#=
function shortest_paths(diags::BitMatrix)::Tuple{Matrix{Int},Matrix{Set{CartesianIndex{2}}}}
    n, m = size(diags)
    d = (0:n) .+ (0:m)'
    b = similar(d, Set{CartesianIndex{2}})
    b[:, 1] .= (eltype(b)(),)
    b[1, :] .= (b[1, 1],)

    for i = 1:n, j = 1:m
        shortest = min(d[i+1, j], d[i, j+1])
        if diags[i, j]
            @assert d[i, j] <= shortest
            if d[i, j] == shortest
                b[i+1, j+1] = intersect(b[i+1, j], b[i, j], b[i, j+1])
            else
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
    d, b
end

lindices(n::Int) = collect(zip([1:n; fill(n, n - 1)], [fill(0, n); 1:n-1]))
tindices(n::Int) = reverse.(lindices(n))

taxicab(v, w)::Int16 = LinearAlgebra.norm(v .- w, 1)
euclid(v, w)::Float64 = LinearAlgebra.norm(v .- w, 2)

#distance_deltas(dg, de) = cumsum(cumsum(@.(2 * (dg - de) - 1), dims = 1), dims = 2)
distance_deltas(dg, de) = @.(2 * (dg - de) - 1)

#=
struct Box
    l::UnitRange{Int}
    t::UnitRange{Int}
end

struct Region
end

function update(r::Region, x::Int, y::Int)
    r.xt < x <= r.xl || return
    r.yl < y <= r.yt || return

    x < r.x

struct Cell
    x::Int
    y::Int
    pre::Box
    boxes::Vector{Box}
    post::Box

    function Box(n::Int, x::Int, y::Int)
        l = x:n+y-1
        t = y:n+x-1
        new(x, y, l, t, [Box(l, t)])
    end
end

function whack()
    if
=#

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
        dd = distance_deltas(dg, de)

        new(n, falses(n, n), dil, dit, dg, de, dd)
    end
end

function add_diag(g::Grid, x::Int, y::Int)::Nothing
    @assert !g.diags[x, y]
    g.diags[x, y] = true

    lix = x:g.n+y-1
    dxyl = g.dil[x, y, lix]
    if x > 1 && y < g.n
        sps, bb = shortest_paths(g.diags[x-1:-1:2, y+1:g.n-1])
        dil = view(g.dil, x-1:-1:1, y+1:g.n, lix)
        dil .= min.(dil, reshape(dxyl, 1, 1, :) .+ sps .+ 1)
    end

    tix = y:g.n+x-1
    dxyt = g.dit[x, y, tix]
    if x < g.n && y > 1
        sps, bb = shortest_paths(g.diags[x+1:g.n-1, y-1:-1:2])
        dit = view(g.dit, x+1:g.n, y-1:-1:1, tix)
        dit .= min.(dit, reshape(dxyt, 1, 1, :) .+ sps .+ 1)
    end

    dg = view(g.dg, lix, tix)
    dg .= min.(dg, reshape(dxyl, :, 1) .+ reshape(dxyt, 1, :) .+ 1)
    g.dd[lix, tix] = distance_deltas(dg, view(g.de, lix, tix))

    return
end

function score(g::Grid, x::Int, y::Int)::Float64
    lix = x:g.n+y-1
    tix = y:g.n+x-1
    dd = view(g.dd, lix, tix)
    dg = view(g.dg, lix, tix)
    sum(dd[reshape(g.dil[x, y, lix], :, 1) .+ reshape(g.dit[x, y, tix], 1, :) .+ 1 .< dg])
end

function score1(g::Grid, x::Int, y::Int)::Float64
    s = 0.0
    for l=x:g.n+y-1, t=y:g.n+x-1
        if g.dil[x, y, l] + g.dit[x, y, t] + 1 < g.dg[l, t]
            s += 1.0
        end
    end
    s
end


function score2(g::Grid, x::Int, y::Int)::Float64
    lix = x:g.n+y-1
    tix = y:g.n+x-1
    sum(reshape(g.dil[x, y, lix], :, 1) .+ reshape(g.dit[x, y, tix], 1, :) .+ 1 .<
        g.dg[lix, tix])
end


function score3(g::Grid, x::Int, y::Int)::Float64
    lix = x:g.n+y-1
    tix = y:g.n+x-1
    dil = view(g.dil, x, y, lix)
    dit = view(g.dit, x, y, tix)
    dg = view(g.dg, lix, tix)
    sum(reshape(dil, :, 1) .+ reshape(dit, 1, :) .+ 1 .< dg)
end

function score4(g::Grid, x::Int, y::Int)::Float64
    lix = x:g.n+y-1
    tix = y:g.n+x-1
    dg = view(g.dg, lix, tix)
    sum(reshape(g.dil[x, y, lix], :, 1) .+ reshape(g.dit[x, y, tix], 1, :) .+ 1 .< dg)
end


function add_best(g::Grid)
    best_score = 0.0, 0
    best_ix = nothing
    for i in CartesianIndices(g.diags)
        g.diags[i] && continue
        x, y = i.I

        if !(x == 1 || y == 1 || x == g.n || y == g.n)
            continue
        end

        s2 = -sum(g.diags[max(1, x-1):min(g.n, x+1),max(1, y-1):min(g.n,y+1)])
        s = score(g, x, y), s2
        if s > best_score
            best_score = s
            best_ix = i
        end
    end
    if best_score[1] > 0
        add_diag(g, best_ix.I...)
    end
    best_score, best_ix
end

using Statistics

all_indices(n::Int) =
    ((x, y) for x in CartesianIndices((0:n, 0:n)), y in CartesianIndices((0:n, 0:n))
         if !(x[1] <= y[1] || x[2] >= y[2]))

function distance(g, v, w)
    x1, y1 = v.I
    x2, y2 = w.I
    diags = g.diags[x1:-1:x2+1,y1+1:y2]
    sps, _ = shortest_paths(diags)
    sps[length(sps)]
end


function scores(g)
    pts = []
    for (i,x) in enumerate(lindices(g.n))
       for (j, y) in enumerate(tindices(g.n))
           if x[1] < y[1] || x[2] > y[2]
               continue
           end
           append!(pts, abs(g.dg[i,j] - g.de[i,j]))
       end
    end
    pts
end

function add_all(g::Grid)
    for i in 1:length(g.diags)
        s, ix = add_best(g)
        if ix == nothing
            return
        end
        pts = scores(g)
        println(s, " ", ix, " ", mean(pts), " ", maximum(pts))
    end
end

function expand(g::Grid)
    big = Grid(g.n * 2)
    for i in CartesianIndices(g.diags)
        if g.diags[i]
            x, y = i.I
            add_diag(big, x * 2 - 1, y * 2 - 1)
            add_diag(big, x * 2, y * 2)
        end
    end
    big
end


function grow(n)
    g = Grid(1)
    add_all(g)
    for i in 1:n-1
        g = expand(g)
        println()
        add_all(g)
    end
    g
end

function expand2(g::Grid)
    big = Grid(g.n + 2)
    for i in CartesianIndices(g.diags)
        if g.diags[i]
            x, y = i.I
            add_diag(big, x + 1, y + 1)
        end
    end
    big
end

function grow2(g, n)
    while g.n < n
        g = expand2(g)
        println("n=",g.n)
        add_all(g)
    end
    g
end

function expand3(g::Grid)
    big = Grid(g.n + 1)
    for i in CartesianIndices(g.diags)
        if g.diags[i]
            x, y = i.I
            add_diag(big, x, y)
        end
    end
    big
end

function grow3(n)
    g = Grid(1)
    add_all(g)
    for i in 1:n-1
        g = expand3(g)
        println("n=",g.n)
        add_all(g)
    end
    g
end

function expand4(g::Grid)
    big = Grid(g.n + 1)
    for i in CartesianIndices(g.diags)
        if g.diags[i]
            x, y = i.I
            add_diag(big, x + 1, y + 1)
        end
    end
    big
end

function grow4(n)
    g = Grid(1)
    add_all(g)
    for i in 1:n-1
        g = expand4(g)
        println("n=",g.n)
        add_all(g)
    end
    g
end


function grow5(n)
    g = Grid(1)
    add_all(g)
    for i in 1:n-1
        if i % 2 == 1
            g = expand3(g)
        else
            g = expand4(g)
        end
        println("n=",g.n)
        add_all(g)
    end
    g
end


#=

function foo(g, x, y)
    @assert !g.diags[x, y]
    lix = x:g.n+y-1
    tix = y:g.n+x-1
    g.dg[lix, tix] .> g.dil[x, y, lix] .+ g.dit[x, y, tix]' .+ 1
end

using Random

function bar(n, m)
    g = Grid(n)
    ix = shuffle(reshape(CartesianIndices(g.diags), :))
    for i in 1:m
        add_diag(g, ix[i].I...)
    end
    foo(g, ix[m+1].I...)
end

function baz(n,m,k)
    for _ in 1:k
        d = bar(n,m)
        for i in CartesianIndices(d)
            x0, y0 = i.I
            for x1 in x0+1:size(d)[1]
                for y1 in y0+1:size(d)[2]
                    if d[x0,y0] == d[x0, y1] == d[x1, y0] == d[x1, y1] == true
                        sub = d[x0:x1,y0:y1]
                        @assert (sum(sub) == length(sub) * d[x0,y0]) d
                    end
                end
            end
        end
    end
end

function baz2(n,m,k)
    for _ in 1:k
        d = bar(n,m)
        x, y = size(d)
        for i in 1:1
            if d[i, 1] == d[i, y] == false
                @assert (sum(d[i, :]) == 0) d
            end
        end
    end
end

function baz3(n,m,k)
    for _ in 1:k
        d = bar(n,m)
        x, y = size(d)
        if d[:,1] == d[:,y] == d[1,:] == d[x,:] == false
            @assert (sum(d) == 0) d
        end
    end
end




0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0
0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0

0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0]

0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0

0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0


function score(g::Grid, i::Impact)::Float64
    s = g.dd[i.l.stop, i.t.stop]
    if i.l.start > 1
        s -= g.dd[i.l.start-1, i.t.stop]
    end
    if i.t.start > 1
        s -= g.dd[i.l.stop, i.t.start-1]
    end
    if i.l.start > 1 && i.t.start > 1
        s += g.dd[i.l.start-1, i.t.start-1]
    end
    s
end

function update_box(g::Grid, b::Box, x::Int, y::Int)

isempty(i::Impact)::Bool = isempty(i.l) || isempty(i.t)

intersects(i::Impact, j::Impact)::Bool =
    !isempty(intersect(i.l, j.l)) && !isempty(intersect(i.t, j.t))


    for i in CartesianIndices(g.diags)
        g.diags[i] && continue
        g.boxes[i] = vcat((update_box(g, b) for b in g.boxes[i])...)
    end

struct Flip
    x::Int
    y::Int
    impacts::Vector{Impact}

    Flip(n::Int, x::Int, y::Int) = new(x, y, [Impact(x:n+y-1, y:n+x-1)])
end

score(g::Grid, f::Flip)::Float64 = sum(score(g, i) for i in f.impacts)

improves(g::Grid, f::Flip, l::Int, t::Int)::Bool =
    g.dil[f.x, f.y, l] + 1 + g.dit[f.x, f.y, t] < g.dg[l, t]

struct Flipper
    g::Grid
    flips::Matrix{Flip}

    Flipper(n::Int) = new(Grid(n), [Flip(n, i.I...) for i in CartesianIndices((n, n))])
end

score(f::Flipper)::Matrix{Float64} = f.flips .|> x -> score(f.g, x)

function flip(f::Flipper, x::Int, y::Int)::Nothing

end

function (x, y)


    if a.x == c.x
end
=#



#=


function apply(g::Grid, ::Add, given::Add)::Float64

    Threads.@threads for ij in  CartesianIndices(scores)
        i, j = ij.I
        #delta = view(deltas, i:n_inner+j-1, j:n_inner+i-1)
        #scores[ij] = sum(delta)
        scores[ij] = sum(
            deltas[d] for d in CartesianIndices((i:n_inner+j-1, j:n_inner+i-1)))
    end

=#
=#
end # module

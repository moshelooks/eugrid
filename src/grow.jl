function norm(x, y)
    #ex = exp(-0.8*x)
    #ey = exp(-0.8*y)
    #(x*ex+y*ey)/(ex + ey)
    min(x, y)
end


function choose_diag!(a::Atom, tl, l, t, br, avoid=nothing)::Bool
    n = length(l)
    br[1] = l[1] + 1
    br[n+1] = t[n] + 1
    br = view(br, 2:n)
    br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1

    !isnothing(avoid) && avoid[a] && return false

    s = 0.0

    a2sq = a[2]^2
    for i in 1:a[1]-1
        dsq = i^2 + a2sq
        w = atan(a[2] / (dsq - 0.25)) / norm(i, a[2])
        de = sqrt(dsq)
        s += (abs(de - br[i]) - abs(de - tl[i] - 1)) * w
    end

    a1sq = a[1]^2

    dsq = a1sq + a2sq
    w = (atan(a[1] / (a[2]-0.5)) - atan((a[1] - 0.5) / a[2])) / norm(a[1], a[2])
    de = sqrt(dsq)
    s += (abs(de - br[a[1]]) - abs(de - tl[a[1]] - 1)) * w

    for i in a[1]+1:a[1]+a[2]-1
        dsq = a1sq + (a[1]+a[2]-i)^2
        w = atan(a[1] / (dsq - 0.25)) / norm(a[1], a[1]+a[2]-i)
        de = sqrt(dsq)
        s += (abs(de - br[i]) - abs(de - tl[i] - 1)) * w
    end

    if s >= 0
        br .= view(tl, 1:n-1) .+ 1
        true
    else
        false
    end
end

function choose_children!(diags, a::Atom, grandparents, parents, children, avoid=nothing)
    Threads.@threads for j in 1:size(children)[2]
        tl = view(grandparents, :, j)
        l = view(parents, :, j)
        t = view(parents, :, j+1)
        br = view(children, :, j)
        a_j = a + CartesianIndex(-1, 1) * (j - 1)
        diags[a_j] = choose_diag!(a_j, tl, l, t, br, avoid)
    end
end

function children!(diags, a::Atom, grandparents, parents, children)
    Threads.@threads for j in 1:size(children)[2]
        tl = view(grandparents, :, j)
        l = view(parents, :, j)
        t = view(parents, :, j+1)
        br = view(children, :, j)
        a_j = a + CartesianIndex(-1, 1) * (j - 1)

        n = length(l)
        br[1] = l[1] + 1
        br[n+1] = t[n] + 1
        br = view(br, 2:n)
        if diags[a_j]
            br .= view(tl, 1:n-1) .+ 1
        else
            br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1
        end
    end
end

diags_buffer(n::Int) = Array{Int, 3}(undef, 2*n+1, n+2, 3)

function grow_base!(buffer)
    n = size(buffer, 2) - 2
    diags = BitMatrix(undef, n, n)

    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        choose_children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end

    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)
    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end

    diags
end

function grow_lower!(buffer, diags, avoid=nothing)
    n = checksquare(diags)
    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)

    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children, avoid)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end
end

function grow_diags(n::Int)::BitMatrix
    buffer = diags_buffer(n+1)
    base = grow_base!(buffer)
    diags = base[(n+1)*onexy:-onexy:2*onexy]
    grow_lower!(buffer, diags)
    diags
end

struct Grid
    diags::BitMatrix
    antidiags::BitMatrix
end

Base.getproperty(g::Grid, s::Symbol) = s === :n ? size(g.diags, 1) : getfield(g, s)

#=
function grow_grid(n::Int, fringe=0)::Grid
    diags = grow_diags(n + 2*fringe)
    antidiags = grow_diags(n + 2*fringe, view(diags, size(diags, 1):-1:1, :))
    Grid(diags, antidiags)
end
=#

#=
function grow_grid(n::Int)::Grid
    buffer = diags_buffer(n+1)
    base = grow_base!(buffer)
    diags = base[(n+1)*onexy:-onexy:2*onexy]
    grow_lower!(buffer, diags)

    avoid = BitMatrix(undef, n+1, n+1)
    avoid[2:n+1, 2:n+1] .= view(diags, :, n:-1:1)
    grow_lower!(buffer, base, avoid)
    antidiags = base[(n+1)*onexy:-onexy:2*onexy]
    grow_lower!(buffer, antidiags, diags[n:-1:1, :])

    Grid(diags, antidiags)
end
=#

function grow_grid(n::Int)::Grid
    buffer = diags_buffer(2n)
    base = grow_base!(buffer)
    diags = base[(n+1)*onexy:2n*onexy]

    avoid = falses(2n, 2n)
    avoid[(n+1)*onexy:2n*onexy] .= view(diags, n:-1:1, :)
    grow_lower!(buffer, base, avoid)
    antidiags = base[(n+1)*onexy:2n*onexy]

    Grid(diags, antidiags)
end

#=
    base = grow_diags!(buffer)
    diags = grow_lower!(buffer, base)

    diags[2*onexy:(n+1)*onexy]
    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        choose_children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)

    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end

    diags = diags[Atom(size(diags)):-onexy:onexy]

    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)

    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children, avoid)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end

    diags
end
=#


function score(g)
    g = BitMatrix(g)
    n = size(g)[1]
    step = Int(n / 8)
    sum(arc(g, step) .!= earc(n, step)) / n^2
end

function grow_simple(n::Int)::BitMatrix
    ds = Matrix{Int}(undef, n+1, n+1)
    ds[:, 1] = ds[1, :] = 0:n
    diags = BitMatrix(undef, n, n)
    for i in CartesianIndices(diags)
        tl = ds[i]
        l = ds[i+onex]
        t = ds[i+oney]
        de = sqrt(i[1]^2 + i[2]^2)
        d_diag = tl + 1
        d_nodiag = min(l, t) + 1
        if abs(de - d_diag) < abs(de - d_nodiag)
            ds[i + onexy] = d_diag
            diags[i] = true
        else
            ds[i + onexy] = d_nodiag
            diags[i] = false
        end
    end
    diags
end

function minscore(g, n, skip=1)
    to = Atom(n-1,n-1)
    ix = onexy:Atom(skip, skip):Atom(size(g))-to
    scores = [score(g[i:i+to]) for i in ix]
    s, i = findmin(scores)
    s, ix[i]
end


#=
struct Geodesic
    lambda::Int
    nodes::Set{Atom}
end


extend(g::Geodesic, a::Atom) = Geodesic(g.lambda + 1, union(g.nodes, [a]))

function nodiag(l::Geodesic, t::Geodesic, a::Atom)
    l.lambda < t.lambda && return extend(l, a)
    t.lambda < l.lambda && return extend(t, a)
    Geodesic(l.lambda + 1, union(l.nodes, t.nodes, [a]))
end

function geodesics(diags)::Matrix{Geodesic}
    ags = Matrix{Geodesic}(undef, size(diags) .+ 1)
    for i in 1:size(diags)[1]+1
        ags[i, 1] = Geodesic(i-1, Set(onexy:Atom(i, 1)))
    end
    for i in 2:size(diags)[2]+1
        ags[1, i] = Geodesic(i-1, Set(onexy:Atom(1, i)))
    end
    for i in onexy:Atom(size(diags))
        br = i + onexy
        ags[br] = diags[i] ? extend(ags[i], br) : nodiag(ags[i+onex], ags[i+oney], br)
    end
    ags
end

function ngeo(
    diags,
    ngeos::Vector{Vector{Int}}=[Vector{Int}() for _ in 1:sum(size(diags))])::
        Vector{Vector{Int}}
    for g in geodesics(diags)
        g.lambda > 0 && push!(ngeos[g.lambda], length(g.nodes))
    end
    ngeos
end

import Statistics

function all_ngeos(diags, k=maximum(sps(diags)))
    kmax = CartesianIndex(k, k)
    ngeos = [Vector{Int}() for _ in 1:sum(kmax.I)]
    for i in CartesianIndices(diags)
        ngeo(view(diags, i:min(i + kmax - onexy, CartesianIndex(size(diags)))), ngeos)
    end
    ngeos[1:k]
end
=#

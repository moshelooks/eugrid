const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

function sps(diags::AbstractMatrix{Bool}, d=Matrix{Int}(undef, size(diags) .+ 1))
    d[:, 1] = 0:size(diags, 1)
    d[1, :] = 0:size(diags, 2)
    for i in CartesianIndices(diags)
        @inbounds d[i + onexy] = 1 + (diags[i] ? d[i] : min(d[i+onex], d[i+oney]))
    end
    d
end

const Distance = NTuple{2, Int32}

b2d(b::Bool)::Distance = (Int32(2) - b, Int32(0))
d2b(d::Distance)::Bool = Int32(2) - d[1]
i2d(i::Integer, j::Integer=0)::Distance = (Int32(i), Int32(j))

d2i(d::Distance, ::RoundingMode{:Down})::Int32 = d[1]
d2i(d::Distance, ::RoundingMode{:Up})::Int32 = d[1] + (d[2] > 0)

disorder(rng::StableRNG, d::Distance)::Distance =
    d[1] == 2 ? d : i2d(1, rand(rng, 1:256))

function sps(dd::AbstractMatrix{Distance}, d=Matrix{Distance}(undef, size(dd) .+ 1))
    d[:, 1] .= i2d.(0:size(dd, 1))
    d[1, :] .= i2d.(0:size(dd, 2))
    for i in CartesianIndices(dd)
        @inbounds d[i + onexy] = min(dd[i] .+ d[i], i2d(1) .+ min(d[i+onex], d[i+oney]))
    end
    d
end

struct Grid
    dd::Matrix{Distance}
    add::Matrix{Distance}

    Grid(::UndefInitializer, n::Int) =
        new(Matrix{Distance}(undef, n-1, n-1), Matrix{Distance}(undef, n-1, n-1))

    function Grid(diags::AbstractMatrix{Bool}, antidiags::AbstractMatrix{Bool})
        checksquare(diags) == checksquare(antidiags) || throw(DimensionMismatch(
            "size(diags) $(size(diags)) != size(antidiags) $(size(antidiags))"))
        new(b2d.(diags), b2d.(antidiags))
    end
end

Grid(diags::AbstractMatrix{Bool}, n::Int=checksquare(diags)) =
    Grid(view(diags, onexy:onexy*n), view(diags, n:-1:1, 1:n))

Base.getproperty(g::Grid, s::Symbol) =
    s === :n ? size(g.dd, 1) + 1 :
    s == :diags ? d2b.(g.dd) :
    s == :antidiags ? d2b.(g.add) :
    getfield(g, s)

chessboard(n::Int) = Grid(trues(n-1, n-1), trues(n-1, n-1))
manhattan(n::Int) = Grid(falses(n-1, n-1), falses(n-1, n-1))

randgrid(rng::StableRNG, n::Int, p::Float64) =
    Grid(rand(Float64, (n-1, n-1)) .< p, rand(Float64, (n-1, n-1)) .< p)

function disorder(rng::StableRNG, g::Grid)::Grid
    disordered = Grid(undef, g.n)
    disordered.dd .= disorder.(rng, g.dd)
    disordered.add .= disorder.(rng, g.add)
    disordered
end

isplanar(g::Grid)::Bool = !any(zip(g.dd, g.add)) do (d, ad); d2b(d) && d2b(ad); end

const Vertex = CartesianIndex{2}

Base.clamp(v::Vertex, g::Grid) = Vertex(clamp.(v.I, 1, g.n))

vertices(n::Int) = CartesianIndices((n, n))
vertices(g::Grid) = vertices(g.n)

function sps(g::Grid, v::Vertex, m=g.n)::Matrix{Distance}
    d = fill(i2d(typemax(Int32)), (g.n, g.n))

    br = clamp(v+onexy*m, g)
    sps(view(g.dd, v:br-onexy), view(d, v:br))

    tl = clamp(v-onexy*m, g)
    sps(view(g.dd, v-onexy:-onexy:tl), view(d, v:-onexy:tl))

    bl = clamp(v+onex*m-oney*m, g)
    sps(view(g.add, v-oney:onex-oney:bl-onex), view(d, v:onex-oney:bl))

    tr = clamp(v-onex*m+oney*m, g)
    sps(view(g.add, v-onex:oney-onex:tr-oney), view(d, v:oney-onex:tr))

    d
end

distance(g::Grid, u::Vertex, v::Vertex) = sps(g, u, maximum(abs.(u.I .- v.I)))[v]

eccentricity(g::Grid, v::Vertex, dv=sps(g, v))::Int = d2i(maximum(dv), RoundDown)

euclidean_eccentricity(n::Int, v::Vertex)::Float64 =
    maximum([sqrt(sum((u - v).I.^2)) for u in onexy:onexy*(n-1):onexy*n])

expected_euclidean_eccentricity(n::Int)::Float64 =
    Statistics.mean(euclidean_eccentricity(n, v) for v in vertices(n))

geodesics(g::Grid, u::Vertex, v::Vertex, du=sps(g, u), dv=sps(g, v, d2i(du[v], RoundDown)))::
    BitMatrix = [du[i] .+ dv[i] == du[v] for i in vertices(g)]

circle_points(g::Grid, v::Vertex, r::Int, dv=sps(g, v, r)) =
    findall(d->d2i(d, RoundDown) == r, dv)

function midpoints(g::Grid, u::Vertex, v::Vertex, du=sps(g, u),
                   dv=sps(g, v, div(d2i(du[v], RoundDown), 2, RoundUp)))
    duv = d2i(du[v], RoundDown)
    half_duv = div(duv, 2)
    pts = intersect!(circle_points(g, u, half_duv, du),
                     circle_points(g, v, duv - half_duv, dv))
    if duv % 2 == 1
        union!(pts, intersect!(circle_points(g, u, half_duv + 1, du),
                               circle_points(g, v, duv - half_duv - 1, dv)))
    end
    pts
end

function two_circle_points(g::Grid, u::Vertex, v::Vertex, r::Int,
                           du=sps(g, u, r+1), dv=sps(g, v, r+1); strict=true)
    ur = circle_points(g, u, r, du)
    vr = circle_points(g, v, r, dv)
    pts = intersect(ur, vr)
    (!isempty(pts) || strict) && return pts
    union!(
        intersect!(ur, union!(circle_points(g, v, r+1, dv), circle_points(g, v, r-1, dv))),
        intersect!(union!(circle_points(g, u, r+1, du), circle_points(g, u, r-1, du)), vr))
end

function euclidean_arcs(n::Int, reps)::BitMatrix
    r = Int(n / reps)
    plane = BitMatrix(undef, r, r)
    for i in CartesianIndices(plane)
        plane[i] = sqrt(i[1]^2+i[2]^2) < r
    end
    repeat(plane, outer=(reps, reps))
end

function diag_arcs(diags::AbstractMatrix{Bool}, reps::Int)::BitMatrix
    r = Int(checksquare(diags) / reps)
    plane = BitMatrix(undef, size(diags))
    for i in onexy:onexy * r:Vertex(size(diags))
        box = i:i+onexy*(r-1)
        plane[box] .= view(sps(view(diags, box)), 2:r+1, 2:r+1) .< r
    end
    plane
end

score_arcs(diags::AbstractMatrix{Bool}, reps::Int=8)::Float64 =
    sum(diag_arcs(diags, reps) .!= euclidean_arcs(checksquare(diags), reps)) / length(diags)

function crisscross(g::Grid, m::Int)::BitMatrix
    pairs = Set{Set{Vertex}}()
    for i in 0:m
        x = 1+div(i*(g.n-1), m)
        for u in Vertex.([(1, x), (x, 1), (g.n, x), (x, g.n)])
            for j in 0:m
                y = 1+div(j*(g.n-1), m)
                for v in Vertex.([(1, y), (y, 1), (g.n, y), (y, g.n)])
                    u != v && push!(pairs, Set{Vertex}([u, v]))
                end
            end
        end
    end
    plane = falses(g.n, g.n)
    for (u, v) in pairs
        plane .|= geodesics(g, u, v)
    end
    plane
end

function mangle(g::Grid, k::Int, v::Vertex, dv=sps(g, v, k))
    dv = d2i.(dv, RoundDown)
    function along((xhat, yhat))
        extent = v + k*xhat
        clamp(extent, g) != extent && return nothing
        clamp(extent + k*yhat, g) != extent + k*yhat && return nothing
        clamp(extent - k*yhat, g) != extent - k*yhat && return nothing
        @assert dv[extent] == k "$(dv[extent]) vs. $k"
        sum(dv[extent+yhat:extent+k*yhat] .== k), sum(dv[extent-k*yhat:extent-yhat] .== k)
    end

    ys = filter(!isnothing, along.([
        (onex, oney), (oney, onex), (-onex, oney), (-oney, onex)]))
    Statistics.mean([atand(a, k) + atand(b, k) for (a, b) in ys])
end

function rmangle(rng::StableRNG, g::Grid, k::Int, m::Int)
    vs = vertices(g)[k+1:g.n-k,k+1:g.n-k]
    [mangle(g, k, rand(rng, vs)) for _ in 1:m]
end

function rmangle(rng::StableRNG, g::Grid, ks::Vector{Int}, m::Int)
    k = maximum(ks)
    vs = vertices(g)[k+1:g.n-k,k+1:g.n-k]
    ys = Matrix{Float64}(undef, m, length(ks))
    for i in 1:m
        v = rand(rng, vs)
        dv = sps(g, v, k)
        for j in eachindex(ks)
            ys[i, j] = mangle(g, ks[j], v, dv)
        end
    end
    ys
end

function explode(g::Grid)
    d2 = falses(2*(g.n-1), 2*(g.n-1))
    ad2 = falses(2*(g.n-1), 2*(g.n-1))
    diags, antidiags = g.diags, g.antidiags
    for i in CartesianIndices(diags)
        if diags[i]
            d2[2*i] = d2[2*i-onexy] = true
        end
        if antidiags[i]
            ad2[2*i-onex] = ad2[2*i-oney] = true
        end
    end
    Grid(d2, ad2)
end

#=
function count_along(xs, start)
    x = xs[start]
    n = 1
    for i in start+1:length(xs)
        xs[i] != x && break
        n += 1
    end
    for i in start-1:-1:1
        xs[i] != x && break
        n += 1
    end
    n
end

function mangle(g::Grid, v::Vertex, dv=sps(g, v))
    dmin = 4
    dmax = round(euclidean_eccentricity(g.n, v), RoundDown)
    euclidean = [Set{Vertex}() for _ in dmin:dmax]
    grid = [Set{Vertex}() for _ in dmin:dmax]
    for i in CartesianIndices(dv)
        ed = isqrt(sum((v.I .- i.I).^2))
        ed < dmin && continue
        push!(euclidean[ed - dmin + 1], i)
        gd = d2i(dv[i], RoundDown)
        (gd < dmin || gd > dmax)  && continue
        push!(grid[gd - dmin + 1], i)
    end
    nums = [length(intersect(e, g)) for (e, g) in zip(euclidean, grid)]
    denoms = [length(union(e, g)) for (e, g) in zip(euclidean, grid)]
    xs = Int[]
    ys = Float64[]
    mind = floor(2pi * dmin)
    for (i, (n, d)) in enumerate(zip(nums, denoms))
        d < mind && continue
        push!(xs, i + dmin - 1)
        push!(ys, n / d)
    end
    #GLM.coef(GLM.lm(@GLM.formula(Y ~ X), DataFrames.DataFrame(X=xs, Y=ys)))[2]
    xs, ys
end

function gruntle(g::Grid, v::Vertex, dv=sps(g, v))
    xs = Float64[]
    ys = Int64[]
    for i in CartesianIndices(dv)
        push!(xs, sqrt(sum((v.I .- i.I).^2)))
        push!(ys, d2i(dv[i], RoundDown))
    end
    xs, ys
end


function mangle_stats(g::Grid, a::Vertex)
    m = Int((g.n - 1) / 2) - 1
    d = [Vector{Int}() for _ in 1:4]
    for i in findall(d2i.(sps(g, a, m), RoundDown) .== m)
        i -= a
        if i[1] > m
            i = Vertex(i[1] - g.n, i[2])
        elseif i[1] < -m
            i = Vertex(i[1] + g.n, i[2])
        end
        if i[2] > m
            i = Vertex(i[1], i[2] - g.n)
        elseif i[2] < -m
            i = Vertex(i[1], i[2] + g.n)
        end
        if i[1] == m
            push!(d[1], i[2])
        elseif i[1] == -m
            push!(d[2], i[2])
        end
        if i[2] == m
            push!(d[3], i[1])
        elseif i[2] == -m
            push!(d[4], i[1])
        end
    end
    d
end

function mangle(g::Grid, a::Vertex)
    m = Int((g.n - 1) / 2) - 1
    d = 0
    for (lo, hi) in extrema.(mangle_stats(g, a))
        d = max(d, min(atand((-lo-1)/m) + atand(hi/m), atand(-lo/m) + atand((hi+1)/m)))
    end
    d
end
#=using Statistics

function rmangle(g::Grid, k, seed=1)
    rng = StableRNG(seed)
    xs = [mangle(g, rand(rng, vertices(g.n))) for _ in 1:k]
    minimum(xs), mean(xs), maximum(xs)
end
=#
=#

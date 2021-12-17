const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

function sps(diags::AbstractMatrix{Bool}, d=Matrix{Int}(undef, size(diags) .+ 1))
    d[:, 1] = 0:size(diags, 1)
    d[1, :] = 0:size(diags, 2)
    for i in CartesianIndices(diags)
        @inbounds d[i + onexy] = 1 + (diags[i] ? d[i] : min(d[i+onex], d[i+oney]))
    end
    d
end

@assert Int === Int64
const Distance = Int64

b2d(b::Bool)::Distance = (2 - b) << 32
d2b(d::Distance)::Bool = 2 - (d >> 32)
i2d(i::Int, j::Int=0)::Distance = (i << 32) + j

d2i(d::Distance, ::RoundingMode{:Down})::Int = d >> 32
d2i(d::Distance, ::RoundingMode{:Up})::Int = (d >> 32) + (mod(d, 2^32) > 0)

disorder(rng::StableRNG, d::Distance)::Distance =
    d2i(d,  RoundDown) == 2 ? d : i2d(1, rand(rng, 1:255))

function sps(dd::AbstractMatrix{Distance}, d=Matrix{Distance}(undef, size(dd) .+ 1))
    d[:, 1] .= i2d.(0:size(dd, 1))
    d[1, :] .= i2d.(0:size(dd, 2))
    for i in CartesianIndices(dd)
        @inbounds d[i + onexy] = min(dd[i] + d[i], i2d(1) + min(d[i+onex], d[i+oney]))
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

Grid(diags::AbstractMatrix{Bool}, n::Int=checksquare(diags)+1) =
    Grid(view(diags, onexy:onexy*(n-1)), view(diags, n-1:-1:1, 1:n-1))

Base.getproperty(g::Grid, s::Symbol) =
    s === :n ? size(g.dd, 1) + 1 :
    s == :diags ? d2b.(g.dd) :
    s == :antidiags ? d2b.(g.add) :
    getfield(g, s)

chessboard(n::Int) = Grid(trues(n-1, n-1), trues(n-1, n-1))
manhattan(n::Int) = Grid(falses(n-1, n-1), falses(n-1, n-1))

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

function sps(g::Grid, v::Vertex, m=g.n)
    br = clamp(v+onexy*m, g)
    tl = clamp(v-onexy*m, g)
    bl = clamp(v+onex*m-oney*m, g)
    tr = clamp(v-onex*m+oney*m, g)

    d = OffsetArrays.OffsetArray(fill(i2d(Int(typemax(Int32))), size(tl:br)), tl:br)

    sps(view(g.dd, v:br-onexy), view(d, v:br))
    sps(view(g.dd, v-onexy:-onexy:tl), view(d, v:-onexy:tl))
    sps(view(g.add, v-oney:onex-oney:bl-onex), view(d, v:onex-oney:bl))
    sps(view(g.add, v-onex:oney-onex:tr-oney), view(d, v:oney-onex:tr))

    d
end

distance(g::Grid, u::Vertex, v::Vertex) = sps(g, u, maximum(abs.(u.I .- v.I)))[v]

eccentricity(g::Grid, v::Vertex, dv=sps(g, v))::Int = d2i(maximum(dv), RoundDown)

euclidean_eccentricity(n::Int, v::Vertex)::Float64 =
    maximum([sqrt(sum((u - v).I.^2)) for u in onexy:onexy*(n-1):onexy*n])

HG(g::Grid, v::Vertex, dv=sps(g, v))::Float64 =
    (eccentricity(g, v, dv) - euclidean_eccentricity(g.n, v)) / log2(g.n^2)
# nb. denominator ought to be sqrt(g.n); gets corrected by draw_hg() in figures.jl

geodesics(g::Grid, u::Vertex, v::Vertex,
          du=sps(g, u, sum(abs.(u.I .- v.I))),
          dv=sps(g, v, sum(abs.(u.I .- v.I))))::BitMatrix =
    [du[i] .+ dv[i] == du[v] for i in intersect(CartesianIndices(du), CartesianIndices(dv))]

NG(g::Grid, u::Vertex, v::Vertex, du=sps(g, u), dv=sps(g, v, d2i(du[v], RoundDown)))::Int =
    sum(geodesics(g, u, v, du, dv))

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

TG(g::Grid, d::Int, dmid::Int)::Float64 = (dmid / d - sqrt(3) / 2)

function TG(rng::StableRNG, g::Grid, u::Vertex, v::Vertex, d::Int,
            du=sps(g, u, d+1), dv=sps(g, v, d+1))
    w = rand(rng, two_circle_points(g, u, v, d, du, dv, strict=false))
    @assert w != v

    m = rand(rng, midpoints(g, u, v, du, dv))
    dmid = d2i(distance(g, m, w), RoundDown)

    TG(g, d, dmid)
end

function graphl(g::Grid, k::Int, xhat::Vertex, yhat::Vertex, u::Vertex, du=sps(g, u, k))::Int
    v = u + k*xhat
    for i in 1:k
        d2i(du[v + i*yhat], RoundDown) > k && return i-1
    end
    k
end

function euclideanl(k::Int)::Int
    k2 = k^2
    for i in 1:k
        d = sqrt(k2+i^2)
        @assert abs(d - k) != abs(d - k - 1)
        abs(d - k) >= abs(d - k - 1) && return i-1
    end
    k
end

midpoints(d::Int) = unique!([(div(d, 2, r1), div(d, 2, r2))
                             for (r1, r2) in ((RoundUp, RoundDown), (RoundDown, RoundUp))])

function AG(rng::StableRNG, g::Grid, d::Int, u::Vertex, du=sps(g, u, d+1))::Float64
    v = u + d*rand(rng, (onex, -onex, oney, -oney))
    TG(rng, g, u, v, d, du)
end


function diag_arcs(diags::AbstractMatrix{Bool}, reps::Int, r::Int)::BitMatrix
    step = Int(checksquare(diags) / reps)
    plane = BitMatrix(undef, size(diags))
    for i in onexy:onexy * step:onexy * step * reps
        box = i:i+onexy*(step-1)
        plane[box] .= sps(view(diags, box))[1:step,1:step] .<= r
    end
    plane
end

function grid_circs(g::Grid, reps::Int, r::Int)::BitMatrix
    step = Int((g.n-1) / reps)
    plane = falses(g.n, g.n)
    d = fill(i2d(Int(typemax(Int32))), (g.n, g.n))
    for i in onexy:onexy * step:onexy * step * reps
        box = i:i+onexy*(step-1)
        center = i+onexy*div(step, 2)
        plane[box] .= sps(g, center, r, d)[box] .<= i2d(r)
    end
    plane
end

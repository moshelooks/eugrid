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

disorder(rng::StableRNG, d::Distance)::Distance = d .+ i2d(0, rand(rng, 1:2^16))

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

function diag_arcs(diags::BitMatrix, reps::Int)::BitMatrix
    r = Int(checksquare(diags) / reps)
    plane = BitMatrix(undef, size(diags))
    for i in onexy:onexy * r:Vertex(size(diags))
        box = i:i+onexy*(r-1)
        plane[box] .= view(sps(view(diags, box)), 2:r+1, 2:r+1) .< r
    end
    plane
end

score_arcs(diags::BitMatrix, reps::Int=8)::Float64 =
    sum(diag_arcs(diags, reps) .!= euclidean_arcs(checksquare(diags), reps)) / length(diags)

function crisscross(g::Grid, m)::BitMatrix
    plane = falses(g.n, g.n)
    step = div(g.n - 1, m)
    for i in onexy:onexy*step:onexy*g.n
        for j in onexy:onexy*step:onexy*g.n
            plane .|= geodesics(g, i, j)
        end
    end
    plane
end

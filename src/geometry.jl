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
#=
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

#=
n_geodesics(g::Grid, u::Vertex, v::Vertex, du=sps(g, u), dv=sps(g, u, du[v]))::Int =


function geodesics(g::Grid, u::Vertex, v::Vertex)::BitMatrix



function nx(g, p)
    g = deepcopy(g)
    for x in (g.diags, g.antidiags), i in CartesianIndices(x)
        if rand() < p
            x[i] = false
        end
    end
    g
end


function step!(s::State, (grandparents, parents, children); W=nothing, sparsity=nothing)
    scores = score!(s, grandparents, parents, children, W=W)
    #sd = Statistics.std(scores)
    #scores .+= randn(length(scores)) * sd * 1e-1
    #if s.position % 2 == 0
    #    sparsity *= 0.85
    #end
    #sparsity ^= mod1(s.position, 3)
    #cutoff = sparsity
    #cutoff = max(sparsity_cutoff(scores, sparsity), 0.0)
    cutoff = isnothing(sparsity) ? 0.0 : max(0.0, sparsity_cutoff(scores, sparsity))
    #=#cutoff = 0.0#min(cutoff, 2.0)
    #if s.position % 2 == 1
    #    cutoff /= 2
    #end
    =#
    #scores .+= abs.(randn(length(scores)) * 1e-4)
    #cutoff -= abs(randn() * sd)

    diag_indices = findall(scores .>= cutoff)
    #diag_indices = findall((scores .>= cutoff) .& (rand(length(scores)) .>= sparsity))
        #min(0.5, 4 / length(scores))))
    s.diags[s.vertices[diag_indices]] .= true

    depth = size(parents, 1) - 1
    Threads.@threads for i in diag_indices
        children[2:depth+1, i] .= view(grandparents, 1:depth, i) .+ 1
    end
end




Base.size(g::Grid) = size(g.diags)
Base.size(g::Grid, dim) = size(g.diags, dim)
Base.getproperty(g::Grid, s::Symbol) = s === :n ? size(g, 1) : getfield(g, s)



Atoms(n::Int) = CartesianIndices((n, n))
Atoms(a::Atom) = CartesianIndices(a.I)
Atoms(g::Grid) = CartesianIndices(g.diags)





function diag(dtl::DistanceMatrix)::DistanceMatrix
    dbr = DistanceMatrix(dtl.a + onexy)
    dbr.data[Atoms(dtl)] .= dtl.data .+ 1
    dbr
end

function nodiag(dt::DistanceMatrix, dl::DistanceMatrix)::DistanceMatrix
    @assert dt.a + onex == dl.a + oney "top $(dt.a) incompatible with left $(dl.a)"
    dbr = DistanceMatrix(dt.a + onex)
    ix = Atoms(dt.a - oney)
    dbr.data[ix] .= min.(view(dt.data, ix), view(dl.data, ix)) .+ 1
    dbr
end

struct Triple
    a::Int
    b::Int
    c::Int
    Triple(a::Int, b::Int) = new(a, b, sqrt(a^2 + b^2))
end

Triple(v::Atom) = Triple(v.I...)

Base.getproperty(t::Triple, s::Symbol) = s === :v ? Atom(t.a, t.b) : getfield(t, s)

all_triples(n::Int, m::Int = n)::Vector{Triple} =
    [Triple(v) for v in Atoms(n, m) if isinteger(sqrt(sum(v.I.^2)))]

struct Box
    br::Atom
    triples::Vector{Triple}
end

Box(n::Int, m::Int, ts::Vector{Triple} = all_triples(n, m)) =
    Box(Atom(n, m) + onexy, ts)

Box(n::Int, ts::Vector{Triple} = all_triples(n)) = Box(n, n, ts)

struct Region
    t::Triple
    u::Atom
end

Base.getproperty(r::Region, s::Symbol) = s === :w ? r.u+r.t.v : getfield(r, s)

Base.in(a::Atom, r::Region)::Bool = all(r.u.I .<= a.I .< r.w.I)

delta_min(r::Region, a::Atom)::Int = maximum((r.w -a).I)
delta_max(r::Region, a::Atom)::Int = sum((r.w - a).I)

delta_distance(r::Region, d::DistanceMatrix)::Int = r.t.c - d.data[r.u]

regions_of_interest(br::Atom, t::Triple, a::Atom) =
    (Region(t, u) for u in max(onexy, a - t.v + onexy):min(a, br - t.v))
regions_of_interest(b::Box, a::Atom) =
    Iterators.flatten(regions_of_interest(b.br, t, a) for t in b.triples)

struct Constraints
    domain::Set{Atom}
    region_clauses::Dict{Region, Vector{Atom}}

    Constraints() = new(Set{Atom}(), Dict{Region, Vector{Atom}}())
end

const _free = [Atom(0,0)]

_isfree(c)::Bool = c === _free
isfree(cs::Constraints, r::Region)::Bool = _isfree(get(cs.region_clauses, r, nothing))

free!(cs::Constraints, r::Region) = setindex!(cs.region_clauses, _free, r)

cover!(cs::Constraints, r::Region) = get!(Vector{Atom}, cs.region_clauses, r)

issatisfiable(cs::Constraints)::Bool = all(!isempty, values(cs.region_clauses))

violations(cs::Constraints)::Vector{Region} =
    [r for (r, c) in cs.region_clauses if isempty(c)]

const Blocker = Dict{Atom, Bool}

#blockers(cs::Constraints, ds, diags::Set{Atom})::Vector{Blocker} =
#    [Blocker(a=>in(a, diags) for a in keys(ds) if a in r) for r in violations(cs)]

#blockers(cs::Constraints, ds, diags::Set{Atom})::Vector{Blocker} =
#    [Blocker(a=>in(a, diags) for a in keys(ds) if a in r || (a + onex) in r || (a + oney in r)) for r in violations(cs)]

blockers(cs::Constraints, ds, diags::Set{Atom})::Vector{Blocker} =
    [Blocker(a=>in(a, diags) for a in keys(ds) if a[1] < r.w[1] && a[2] < r.w[2])
     for r in violations(cs)]


clauses(cs::Constraints) = Iterators.filter(!_isfree, values(cs.region_clauses))

lowermost(r::Region, a::Atom) = Atom(r.u[1] + r.t.a - 1, a[2])
rightmost(r::Region, a::Atom) = Atom(a[1], r.u[2] + r.t.b - 1)

function constrain(cs::Constraints, r::Region, d::DistanceMatrix, ribbon::AbstractSet{Atom})
    delta = delta_distance(r, d)
    delta == delta_max(r, d.a) && return :negates
    isfree(cs, r) && return nothing
    dmin = delta_min(r, d.a)
    delta > dmin && return :frees
    for (i, x) in enumerate(d.a+onex:lowermost(r, d.a))
        !in(x, ribbon) && delta >= i + delta_min(r, x) && return :frees
    end
    for (i, y) in enumerate(d.a+oney:rightmost(r, d.a))
        !in(y, ribbon) && delta >= i + delta_min(r, y) && return :frees
    end
    delta == dmin && return :satisfies
    return :unsatisfiable
end

function constrain!(cs::Constraints, b::Box, d::DistanceMatrix, ribbon::AbstractSet{Atom})
    negated = false
    clausal_regions = Vector{Region}()
    for r in regions_of_interest(b, d.a)
        constraint = constrain(cs, r, d, ribbon)
        isnothing(constraint) && continue
        if constraint === :negates
            negated = true
            free!(cs, r)
        elseif constraint == :frees
            free!(cs, r)
        elseif constraint == :satisfies && !negated
            push!(clausal_regions, r)
        else
            cover!(cs, r)
        end
    end
    if negated
        foreach(clausal_regions) do r; cover!(cs, r) end
    else
        push!(cs.domain, d.a)
        foreach(clausal_regions) do r; push!(cover!(cs, r), d.a) end
    end
end

const GraphDistances = Dict{Atom, DistanceMatrix}

function constraints(b::Box, ds::GraphDistances)::Constraints
    cs = Constraints()
    for d in values(ds)
        constrain!(cs, b, d, keys(ds))
    end
    filter!(cs.region_clauses) do (r, c)
        c == _free && return false
        if ((haskey(ds, r.w - onex) && delta_distance(r, ds[r.w - onex]) == 1) ||
            (haskey(ds, r.w - oney) && delta_distance(r, ds[r.w - oney]) == 1))
            false
        else
            true
        end
    end
    cs
end

const Ribbon = Vector{Atom}

function ribbons(kernel::Ribbon, basis::Ribbon, n::Int)::Vector{Ribbon}
    ribbons = Vector{Ribbon}()
    atoms = Set(Atoms(n))
    while !isempty(atoms)
        kernel = sort!(filter(a->!isnothing(pop!(atoms, a, nothing)), kernel))
        @assert !isempty(kernel) "incomplete kernel/basis does not cover $atoms"
        push!(ribbons, kernel)
        kernel = map(sum, Iterators.product(kernel, basis))
    end
    ribbons
end

struct Membrane
    interior::GraphDistances
    border::GraphDistances
    targets::Ribbon
    diag_cache::Dict{Atom, DistanceMatrix}
    nodiag_cache::Dict{NTuple{2, UInt64}, DistanceMatrix}

    function Membrane(interior::GraphDistances, exterior::Ribbon)
        border = GraphDistances()
        targets = filter(exterior) do a
            minimum(a.I) > 1 && return true
            border[a] = DistanceMatrix(a)
            false
        end
        diag_cache = Dict{Atom, DistanceMatrix}()
        nodiag_cache = Dict{NTuple{2, UInt64}, DistanceMatrix}()
        new(interior, border, targets, diag_cache, nodiag_cache)
    end
end

function exterior_distances(m::Membrane, diags::Set{Atom})::GraphDistances
    ds = copy(m.border)
    for target in m.targets
        tl = target - onexy
        ds[target] = if tl in diags
            get!(m.diag_cache, target) do; diag(m.interior[tl]) end
        else
            t = target - onex
            dt = get(m.interior, t) do; ds[t] end
            l = target - oney
            dl = get(m.interior, l) do; ds[l] end
            get!(m.nodiag_cache, (objectid(dt), objectid(dl))) do; nodiag(dt, dl) end
        end
    end
    ds
end
=#
=#

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

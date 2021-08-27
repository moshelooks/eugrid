const Atom = CartesianIndex{2}

struct DistanceMatrix
    data::Matrix{Int}

    function DistanceMatrix(a::Atom)
        n, m = a.I
        data = Matrix{Int}(undef, n, m)
        data[1:n, m] = n-1:-1:0
        data[n, 1:m-1] = m-1:-1:1
        new(data)
    end
end

Base.size(d::DistanceMatrix) = size(d.data)
Base.getproperty(d::DistanceMatrix, s::Symbol) = s === :a ? Atom(size(d)) : getfield(d, s)

Atoms(x, y) = CartesianIndices((x, y))
Atoms(n::Int) = CartesianIndices((n, n))
Atoms(a::Atom) = CartesianIndices(a.I)
Atoms(d::DistanceMatrix) = CartesianIndices(d.data)

const onex, oney, onexy = Atoms(0:1, 0:1)[2:end]

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

blockers(cs::Constraints, ds, diags::Set{Atom})::Vector{Blocker} =
    [Blocker(a=>in(a, diags) for a in keys(ds) if a in r) for r in violations(cs)]

clauses(cs::Constraints) = Iterators.filter(!_isfree, values(cs.region_clauses))

function constrain!(cs::Constraints, b::Box, d::DistanceMatrix, wiggle=0)::Nothing
    negated = false
    clausal_regions = Vector{Region}()
    for r in regions_of_interest(b, d.a)
        delta = delta_distance(r, d)
        if delta == delta_max(r, d.a) + wiggle
            negated = true
            free!(cs, r)
        elseif !isfree(cs, r)
            delta -= delta_min(r, d.a)
            if delta > -wiggle
                free!(cs, r)
            elseif delta == -wiggle && !negated
                push!(clausal_regions, r)
            else
                cover!(cs, r)
            end
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
        constrain!(cs, b, d)
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

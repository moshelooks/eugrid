const DistanceMatrix = Matrix{Int}
const GraphDistances = Vector{DistanceMatrix}

struct Triple
    a::Int
    b::Int
    c::Int
    Triple(a::Int, b::Int) = new(a, b, sqrt(a^2 + b^2))
end

vrange(t::Triple, k::Int, m::Int) = max(k, m-t.a+1):m-1
hrange(t::Triple, k::Int, m::Int) = max(m-t.a >= k ? 1 : k, m-t.b+1):m-1

RegionIndices(t::Triple, k::Int, m::Int) =
    unique(Iterators.flatten((CartesianIndices((vrange(t, k, m), k:k)),
                              CartesianIndices((k:k, hrange(t, k, m))))))

struct Region
    t::Triple
    v::UnitRange{Int}
    h::UnitRange{Int}

    function Region(t::Triple, u::CartesianIndex{2}, m::Int)
        vstart, hstart = u.I
        @assert max(vstart, hstart) <= m
        vstop, hstop = min.((vstart, hstart) .+ (t.a, t.b) .- 1, m)
        @assert max(vstop, hstop) >= m
        if vstop < m
            hstop = hstart-1
        elseif hstop < m
            vstop = vstart-1
        end
        new(t, vstart:vstop, hstart:hstop)
    end
end

vorigin(r::Region)::CartesianIndex{2} = CartesianIndex(r.v.start, r.h.start)
horigin(r::Region)::CartesianIndex{2} = CartesianIndex(r.h.start, r.v.start)

lastrow(r::Region)::Int = r.v.start + r.t.a
lastcol(r::Region)::Int = r.h.start + r.t.b

lastk(r::Region)::Int =
    isempty(r.v) ? r.h.stop : isempty(r.h) ? r.v.stop : max(r.v.stop, r.h.stop)

regions(t::Triple, k::Int, m::Int) = (Region(t, u, m) for u in RegionIndices(t, k, m))

negation_distance(r::Region, k::Int, m::Int)::Int = r.t.c + k + m - lastrow(r) - lastcol(r)

function negates(r::Region, d::DistanceMatrix)::Bool
    k, m = size(d)
    target = negation_distance(r, k, m)
    k in r.v && d[vorigin(r)] == target && return true
    k in r.h && d[horigin(r)] == target && return true
    false
end

affirmation_distance(r::Region, k::Int, m::Int)::Int =
    r.t.c - max(lastrow(r) - k, lastcol(r) - m)

function affirmation_bound(r::Region, d::DistanceMatrix)::Int
    k, m = size(d)
    if k in r.v
        bound = d[vorigin(r)] - affirmation_distance(r, k, m)
        if bound >= 0 && k in r.h
            bound = min(bound, d[horigin(r)] - affirmation_distance(r, m, k))
        end
    else
        @assert k in r.h "$k $m $r"
        bound = d[horigin(r)] - affirmation_distance(r, m, k)
    end
    bound
end

mutable struct Constraint
    region::Region
    clause::Union{Vector{Int}, Nothing}
end

function Constraint(r::Region, d::DistanceMatrix, negated::Bool)
    bound = affirmation_bound(r, d)
    bound < 0 && return Constraint(r, nothing)
    c = Vector{Int}()
    bound == 0 && !negated && push!(c, size(d)[1])
    Constraint(r, c)
end

function update_clause!(c::Constraint, d::DistanceMatrix, negated::Nothing)::Bool
    !negates(c.region, d) && return false
    c.clause = nothing
    true
end

function update_clause!(c::Constraint, d::DistanceMatrix, negated::Bool)::Nothing
    isnothing(c.clause) && return
    bound = affirmation_bound(c.region, d)
    if  bound < 0
        c.clause = nothing
    elseif bound == 0 && !negated
        push!(c.clause, size(d)[1])
    end
    nothing
end

constraints(ts::Vector{Triple}, d::DistanceMatrix, negated::Bool) =
    (Constraint(r, d, negated) for t in ts for r in regions(t, size(d)...))

struct MonoCNF
    clauses::Vector{Set{Int}}
end

struct UnsatisfiableException <: Exception
    c::Constraint
end

function MonoCNF(ts::Vector{Triple}, ds::GraphDistances)::MonoCNF
    clauses = Vector{Vector{Int}}()
    active_constraints = Vector{Constraint}()
    for (k, d) in enumerate(ds)
        negated = mapreduce(|, active_constraints, init=false) do c
            update_clause!(c, d, nothing)
        end
        filter!(active_constraints) do c
            update_clause!(c, d, negated)
            k != lastk(c.region) && return true
            if !isnothing(c.clause)
                isempty(c.clause) && throw(UnsatisfiableException(c))
                push!(clauses, c.clause)
            end
            false
        end
        append!(active_constraints, constraints(ts, d, negated))
    end
    MonoCNF(Set{Int}.(sort!(clauses)))
end

literal_counts(cnf::MonoCNF)::Dict{Int, Int} =
    mapreduce(c->Dict(c.=>1), mergewith(+), cnf.clauses)

function simplify!(cnf::MonoCNF)::MonoCNF
    for subset in sort(cnf.clauses, by=length)
        filter!(cnf.clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end
    cnf
end

function solve!(cnf::MonoCNF)::MonoCNF
    isempty(cnf.clauses) && return cnf
    unique!(cnf.clauses)
    while true
        simplify!(cnf)
        n, l = maximum(reverse, literal_counts(cnf))
        n == 1 && break
        push!(cnf.clauses, Set((l,)))
    end
    for c in cnf.clauses
        filter!(>=(maximum(c)), c)
    end
    sort!(cnf.clauses, by=only)
    cnf
end

function solution(ts::Vector{Triple}, ds::GraphDistances)::BitVector
    diags = falses(length(ds))
    diags[map(only, solve!(MonoCNF(ts, ds)).clauses)] .= true
    diags
end

function graph_distances(dps::GraphDistances, diags::AbstractVector{Bool})::GraphDistances
    n = length(dps)
    @assert length(diags) == n
    ds = DistanceMatrix.(undef, 1:n+1, n+1)
    ds[1][1,:] .= n:-1:0
    for (k, diag) in enumerate(diags)
        @views @. begin
            ds[k+1][1:k,1:n] = 1 + (
                diag ? dps[k] : k < n ? min(ds[k][:, 1:n], dps[k+1][1:k,:]) : ds[k][:, 1:n])
            ds[k+1][k+1,1:n] = n:-1:1
            ds[k+1][:,n+1] = k:-1:0
        end
    end
    ds
end

function grow(ts::Vector{Triple}, n::Int)::BitMatrix
    eugrid = falses(n, n)
    ds = GraphDistances()
    diags = BitVector()
    for m in 1:n
        ds = graph_distances(ds, diags)
        eugrid[1:m, m] = eugrid[m, 1:m] = diags = solution(ts, ds)
    end
    eugrid
end

function grow(n::Int, m::Int=n)::BitMatrix
    ts = Vector{Triple}()
    for b in 1:n
        for a in 1:b
            isinteger(sqrt(a^2+b^2)) && push!(ts, Triple(a,b))
        end
    end
    println(ts)
    grow(ts, m)
end

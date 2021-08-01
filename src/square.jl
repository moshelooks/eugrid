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


struct UnsatisfiableException <: Exception
    c::Constraint
end

struct BacktrackingException <: Exception
end

include("monocnf.jl")


function MonoCNF(ts::Vector{Triple}, ds::GraphDistances)::MonoCNF
    clauses = Vector{Vector{Int}}()
    free = Set{Int}()
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
        !negated && push!(free, k)
    end
    MonoCNF(Set{Int}.(sort!(clauses)), free, Vector{Int}())
end

literal_counts(cnf::MonoCNF)::Dict{Int, Int} =
    mapreduce(c->Dict(c.=>1), mergewith(+), cnf.clauses)

function simplify!(cnf::MonoCNF)::MonoCNF
    unique!(cnf.clauses)
    for subset in sort(cnf.clauses, by=length)
        filter!(cnf.clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end
    cnf
end

function solve!(cnf::MonoCNF)::MonoCNF
    isempty(cnf.clauses) && return cnf
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

function affirm_singletons!(cnf::MonoCNF)::MonoCNF
    filter!(cnf.clauses) do c
        length(c) > 1 && return true
        l = only(c)
        pop!(cnf.free, l)
        push!(cnf.affirmed, l)
        false
    end
    cnf
end

function choose!(cnf::MonoCNF)::MonoCNF
    println("choose $cnf")

    l = pop!(cnf.free)

    affirmed_clauses = [copy(c) for c in cnf.clauses if !(l in c)]
    affirmed_cnf = MonoCNF(affirmed_clauses, copy(cnf.free), push!(copy(cnf.affirmed), l))

    println("XXX $l $cnf)")
    for c in cnf.clauses
        pop!(c, l, nothing)
    end
    println("YYY $l $cnf")
    validate(cnf)
    affirm_singletons!(cnf)
    validate(cnf)

    affirmed_cnf
end

struct Solver
    stack::Vector{MonoCNF}

    Solver(cnf::MonoCNF) = new([affirm_singletons!(simplify!(deepcopy(cnf)))])
end

Base.isempty(solver::Solver)::Bool = isempty(solver.stack)

function Base.popfirst!(solver::Solver)::Vector{Int}
    top = pop!(solver.stack)
    while !isempty(top.free)
        push!(solver.stack, choose!(top))
    end
    top.affirmed
end

Base.iterate(solver::Solver, state=nothing) =
    isempty(solver) ? nothing : (popfirst!(solver), nothing)

Base.eltype(::Type{Solver}) = Vector{Int}

Base.IteratorSize(::Type{Solver}) = Base.SizeUnknown()

function findsol(solver, m, dp, ts)
    while true
        isempty(solver) && throw(BacktrackingException())
        diags = falses(m)
        diags[popfirst!(solver)] .= true
        ds = graph_distances(dp, diags)
        try
            return (ds, diags, MonoCNF(ts, ds))
        catch e
            !isa(e, UnsatisfiableException) && throw(e)
        end
    end
end

function growbt(ts::Vector{Triple}, n::Int)::BitMatrix
    eugrid = falses(n, n)
    ds = graph_distances(GraphDistances(), BitVector())
    solver = Solver(MonoCNF(ts, ds))
    for m in 1:n-1
        ds, diags, cnf = findsol(solver, m, ds, ts)
        eugrid[1:m, m] = eugrid[m, 1:m] = diags
        dp = ds
        solver = Solver(cnf)
        println(cnf)
    end
    diags = falses(n)
    diags[popfirst!(solver)] .= true
    eugrid[1:n, n] = eugrid[n, 1:n] = diags
    eugrid
end

struct State
    diags::BitVector
    ds::GraphDistances
    solver::Solver
end

function State(ts::Vector{Triple})
    diags = BitVector()
    ds = graph_distances(GraphDistances(), diags)
    solver = Solver(MonoCNF(ts, ds))
    State(diags, ds, solver)
end

function next(ts::Vector{Triple}, s::State)
    ds, diags, cnf = findsol(s.solver, length(s.diags)+1, s.ds, ts)
    State(diags, ds, Solver(cnf))
end

function growbtbt(ts::Vector{Triple}, n::Int)
    stack = [State(ts)]
    while length(stack) < n
        try
            println(stack[end].solver)
            push!(stack, next(ts, stack[end]))
        catch e
            !isa(e, BacktrackingException) && throw(e)
            pop!(stack)
        end
    end
    eugrid = falses(n, n)
    for m in 1:n-1
        eugrid[1:m, m] = eugrid[m, 1:m] = stack[m+1].diags
    end
    diags = falses(n)
    diags[popfirst!(stack[n].solver)] .= true
    eugrid[1:n, n] = eugrid[n, 1:n] = diags
    eugrid
end

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
#    CartesianIndex(k, max(k+t.a-1, m-t.b+1)):CartesianIndex(k, m-1)

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
#horigin(r::Region)::CartesianIndex{2} = CartesianIndex(-1, -1)

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

function pose(ts::Vector{Triple}, ds::GraphDistances)::Union{MonoCNF, Nothing}
    #length(ds) < 3 && return MonoCNF()
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
                push!(clauses, c.clause)
            end
            false
        end
        for c in clauses
            isempty(c) && return nothing
        end
        append!(active_constraints, constraints(ts, d, negated))
        #k == length(ds) && break
        !negated && push!(free, k)
    end
    MonoCNF(Set{Int}.(sort!(clauses)), free, Vector{Int}())
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

function alltriples(n)
    ts = Vector{Triple}()
    for b in 1:n
        for a in 1:b
        #for a in 1:n
            isinteger(sqrt(a^2+b^2)) && push!(ts, Triple(a,b))
        end
    end
    ts
end

function literal_counts(cnf::MonoCNF)::Dict{Int, Int}
    counts = Dict(cnf.free.=>0)
    for c in cnf.clauses
        for l in c
            counts[l] += 1
        end
    end
    counts
end


function choose!(cnf::MonoCNF)
    #println("choose $cnf")

    #n, l = maximum(reverse, literal_counts(cnf))
    #pop!(cnf.free, l)
    l = pop!(cnf.free, maximum(cnf.free))

    affirmed_clauses = [copy(c) for c in cnf.clauses if !(l in c)]
    affirmed_cnf = MonoCNF(affirmed_clauses, copy(cnf.free), push!(copy(cnf.affirmed), l))

    #println("XXX $l $cnf)")
    for c in cnf.clauses
        pop!(c, l, nothing)
    end
    #println("YYY $l $cnf")
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
    #top = popat!(solver.stack, argmin([length(x.free) for x in solver.stack]))
    while !isempty(top.free)
        #push!(solver.stack, choose!(top))
        #push!(solver.stack, top)
        #top = popat!(solver.stack, argmin([length(x.free) for x in solver.stack]))
        tmp = choose!(top)
        push!(solver.stack, top)
        top = tmp
        #tmp, n = choose!(top)
        #if n > 0
        #    push!(solver.stack, tmp)
        #else
        #    push!(solver.stack, top)
        #    top = tmp
        #end
    end
    top.affirmed
end

Base.iterate(solver::Solver, state=nothing) =
    isempty(solver) ? nothing : (popfirst!(solver), nothing)

Base.eltype(::Type{Solver}) = Vector{Int}

Base.IteratorSize(::Type{Solver}) = Base.SizeUnknown()

function findsol(solver, m, dp, ts)
    while true
        isempty(solver) && return nothing
        diags = falses(m)
        diags[popfirst!(solver)] .= true
        ds = graph_distances(dp, diags)
        cnf = pose(ts, ds)
        !isnothing(cnf) && return (ds, diags, cnf)
    end
end

struct State
    diags::BitVector
    ds::GraphDistances
    solver::Solver
end

function State(ts::Vector{Triple})
    diags = BitVector()
    ds = graph_distances(GraphDistances(), diags)
    solver = Solver(pose(ts, ds))
    State(diags, ds, solver)
end

function next(ts::Vector{Triple}, s::State)
    sol = findsol(s.solver, length(s.diags)+1, s.ds, ts)
    isnothing(sol) && return nothing
    ds, diags, cnf = sol
    State(diags, ds, Solver(cnf))
end

function growbtbt(n::Int, ts::Vector{Triple}=alltriples(n))
    stack = [State(ts)]
    while length(stack) < n
        top = next(ts, stack[end])
        if isnothing(top)
            pop!(stack)
        else
            push!(stack, top)
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

function sps(diags)
    d = Matrix{Int}(undef, size(diags) .+ 1)
    d[:, 1] = 0:size(diags)[1]
    d[1, :] = 0:size(diags)[2]
    for i in CartesianIndices(diags)
        d[i + onexy] = 1 + (diags[i] ? d[i] : min(d[i+onex], d[i+oney]))
    end
    d[2:end, 2:end]
end

function check(diags::BitMatrix, ts::Vector{Triple}=alltriples(maximum(size(diags))))
    for i in CartesianIndices(diags)
        d = sps(diags[onexy:i])
        x, y = size(d)
        for t in ts
            if t.a <= x && t.b <= y
                @assert d[t.a, t.b] == t.c
            end
            if t.b <= x && t.a <= y
                @assert d[t.b, t.a] == t.c
            end
        end
    end
end

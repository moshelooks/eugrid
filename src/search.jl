struct Assignment
    clauses::Vector{Set{Atom}}
    free::Set{Atom}
    affirmed::Set{Atom}
end

Assignment(cs::Constraints)::Assignment =
    _simplify!(Assignment(map(Set{Atom}, clauses(cs)), copy(cs.domain), Set{Atom}()))

function _simplify!(a::Assignment)::Assignment
    unique!(a.clauses)
    for subset in sort(a.clauses, by=length)
        filter!(a.clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end

    filter!(a.clauses) do c
        length(c) > 1 && return true
        l = only(c)
        pop!(a.free, l)
        push!(a.affirmed, l)
        false
    end
    a
end

function _fork!(a::Assignment, l::Atom)::Assignment
    forked_clauses = [copy(c) for c in a.clauses if !in(l, c)]
    forked_free = copy(a.free)
    forked_affirmed = push!(copy(a.affirmed), l)

    foreach(c->delete!(c, l), a.clauses)
    _simplify!(a)

    Assignment(forked_clauses, forked_free, forked_affirmed)
end

fork!(a::Assignment)::Assignment = _fork!(a, pop!(a.free))
fork!(a::Assignment, l::Atom)::Assignment = _fork!(a, pop!(a.free, l))

struct Solver
    stack::Vector{Assignment}
end

Solver(cs::Constraints)::Solver = Solver([Assignment(cs)])

Base.isempty(s::Solver)::Bool = isempty(s.stack)

function Base.popfirst!(s::Solver)::Set{Atom}
    top = pop!(s.stack)
    while !isempty(top.free)
        push!(s.stack, fork!(top))
    end
    top.affirmed
end
#=
struct Layer
    diags::Set{Atom}
    ds::GraphDistances
    solver::Solver

    Layer(diags::Set{Atom}, ds::GraphDistances, cs::Constraints) = new(diags, ds, Solver(cs))
end

function _wrap(m::Membrane, diags::Set{Atom}, ts::Vector{Triple})::Union{Layer, Nothing}
    ds = exterior_distances(m, diags)
    cs = constraints(ts, ds)
    issatisfiable(cs) ? Layer(diags, ds, cs) : nothing
end

function Layer(r::Ribbon)::Layer
    m = Membrane(GraphDistances(), r)
    @assert isempty(m.targets)
    _wrap(m, Set{Atom}(), Vector{Triple}())
end

function wrap!(l::Layer, r::Ribbon, ts::Vector{Triple})::Union{Layer, Nothing}
    m = Membrane(l.ds, r)
    while true
        wrapped = _wrap(m, popfirst!(l.solver), ts)
        !isnothing(wrapped) && return wrapped
        isempty(l.solver) && return nothing
        #revise!(l.solver, constraints.violations)
    end
end
=#
struct Onion
    ribbons::Vector{Ribbon}
    membranes::Vector{Membrane}
    solvers::Vector{Solver}
    diags::Vector{Set{Atom}}
    ts::Vector{Triple}

    function Onion(rs::Vector{Ribbon}, ts::Vector{Triple})
        @assert length(rs) > 1 "need at least two ribbons, got $length(rs)"
        m = Membrane(GraphDistances(), rs[1])
        @assert isempty(m.targets)
        ds = exterior_distances(m, Set{Atom}())
        new(rs, [Membrane(ds, rs[2])], [Solver(constraints(ts, ds))], Set{Atom}[], ts)
    end
end


Onion(kernel::Ribbon, basis::Ribbon, n::Int, ts::Vector{Triple} = all_triples(n)) =
    Onion(ribbons(kernel, basis, n), ts)

function step!(o::Onion)
    diags = popfirst!(o.solvers[end])
    ds = exterior_distances(o.membranes[end], diags)
    cs = constraints(o.ts, ds)
    if issatisfiable(cs)
        push!(o.diags, diags)
        push!(o.solvers, Solver(cs))
        if length(o.solvers) == length(o.ribbons)
            push!(o.diags, popfirst!(o.solvers[end]))
            return true
        end
        push!(o.membranes, Membrane(ds, o.ribbons[length(o.solvers) + 1]))
    else
        while isempty(o.solvers[end])
            pop!(o.membranes)
            pop!(o.solvers)
            isempty(o.diags) && return false
            pop!(o.diags)
        end
    end
    nothing
end

solve!(o::Onion) =
    while true
        s = step!(o)
        !isnothing(s) && return s
    end

function eugrid(o::Onion)
    diags = Matrix{Union{Missing, Bool}}(missing, maximum(o.ribbons[end]).I)
    for (r, d) in zip(o.ribbons, o.diags)
        diags[r] .= in.(r, Ref(d))
    end
    diags
end

#fixme need first ribbon for building eugrid


#=
    m = o.membranes[end]
    s

    next_layer = wrap!(o.layers[end], next_ribbon(o), o.ts)
    if isnothing(next_layer)
        while true




Base.isempty(o::Onion) = isempty(o.layers)


function wrap!(o::Onion)
    n = length(o.layers)
    while length(o.layers) <= n
        l = wrap!(o.layers[end], o.ribbons[length(o.layers)], o.ts)
        if isnothing(l)
            fixme need to keep membrane in the layer
        push!(o.layers, l)

    while





        function eugrid(o::Onion)

Base.size(o::Onion)::NTuple{2, Int} = maximum(o.ribbons[end-1]).I

    diags = Matrix{Union{Missing, Bool}}(missing, size(o))
    for a in Iterators.flatten(l.diags for l in o.layers)
        diags[a] = true
    end
    diags
end

function solve!(o::Onion)
    while !isfull(o)
        l = wrap!(o.layers[end], o.ribbons[length


        wrap!(o) || unwrap
=#

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

function all_solutions!(s::Solver)::Vector{Set{Atom}}
    solutions = Set{Atom}[]
    while !isempty(s)
        push!(solutions, popfirst!(s))
    end
    solutions
end

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

Base.isempty(o::Onion) = all(isempty, o.solvers)

function step!(o::Onion)
    length(o.diags) == length(o.ribbons) && pop!(o.diags)

    while true
        while isempty(o.solvers[end])
            length(o.membranes) == length(o.solvers) && pop!(o.membranes)
            pop!(o.solvers)
            pop!(o.diags)
            isempty(o.solvers) && return false
        end

        diags = popfirst!(o.solvers[end])
        if length(o.solvers) == length(o.ribbons)
            push!(o.diags, diags)
            return true
        end

        ds = exterior_distances(o.membranes[end], diags)
        cs = constraints(o.ts, ds)
        !issatisfiable(cs) && continue
        push!(o.diags, diags)
        push!(o.solvers, Solver(cs))
        length(o.solvers) < length(o.ribbons) &&
            push!(o.membranes, Membrane(ds, o.ribbons[length(o.solvers) + 1]))
        return false
    end
end

function all_steps!(o::Onion)
    steps = Iterators.map(Iterators.countfrom()) do _
        step!(o)
        eugrid(o)
    end
    Iterators.takewhile(_->!isempty(o), steps)
end


function solve!(o::Onion)
    while !isempty(o)
        step!(o) && return true
    end
    false
end

function all_solutions!(o::Onion)
    solutions = BitMatrix[]
    while !isempty(o)
        step!(o) && push!(solutions, eugrid(o))
    end
    solutions
end


function eugrid(o::Onion)
    diags = Matrix{Union{Missing, Bool}}(missing, maximum(o.ribbons[end]).I)
    for (r, d) in zip(o.ribbons, o.diags)
        diags[r] .= in.(r, Ref(d))
    end
    diags
end
#=
function all_solutions(kernel, basis, n)
    o = Onion(kernel, basis, n)
    solutions = BitMatrix[]
    while true
        !solve!(o) && return solutions
        push!(solutions, eugrid(o))
        pop!(o.diags)
        while !isempty(o.solvers[end])
=#

#=

    length(o.diags) == length(o.ribbons) && pop!(o.diags)
    while true
        pop!(o.solvers)
        isempty(o.diags) && return false
        pop!(o.diags)
        length(o.membranes) == length(o.solvers) && pop!(o.membranes)
    end

    else
        pop!(o.membranes)
    end

        pop!(o.diags)
    end




    else
        pop!(o.membran


        pop!(o.membranes)
    while true
        pop!(o.solvers)
        pop!(o.diags)
    if length(o.diags) == length(o.ribbons)


    pop!(o.solvers)


    if !isempty(o.solvers)
        length(o.diags) < length(o.membranes) && pop!(o.membranes)
        pop!(o.diags)
    end
    false
end


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
=#

struct Assignment
    clauses::Vector{Set{Atom}}
    free::Set{Atom}
    affirmed::Set{Atom}
end

Assignment(cs::Constraints)::Assignment =
    _simplify!(Assignment(map(Set{Atom}, clauses(cs)), copy(cs.domain), Set{Atom}()))

function isblocked(a::Assignment, b::Blocker)::Bool
    #=
    for l in a.free
        haskey(b, l) && return false
    end

    for l in a.affirmed
        !get(b, l, true) && return false
    end
    true
    =#
    for (k, v) in b
        k in a.free && return false
        !v && k in a.affirmed && return false
    end
    true
end

isblocked(a::Assignment, bs::Vector{Blocker})::Bool = any(b->isblocked(a, b), bs)

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
    blockers::Vector{Blocker}

    Solver(stack, blockers=Blocker[]) = new(stack, blockers)
end

Solver(cs::Constraints)::Solver = Solver([Assignment(cs)], Blocker[])

state_size(s::Solver) = sum([2^length(a.free) for a in s.stack])

function Base.isempty(s::Solver)::Bool
    isempty(s.stack) && return true
    top = pop!(s.stack)
    while !isempty(top.free)
        l = Atom(minimum(l->l.I, top.free))
        forked = fork!(top, l)
        if isblocked(top, s.blockers)
            if isblocked(forked, s.blockers)
                isempty(s.stack) && return true
                top = pop!(s.stack)
            else
                top = forked
            end
        elseif !isblocked(forked, s.blockers)
            push!(s.stack, forked)
        end
    end
    push!(s.stack, top)
    false
end

function Base.popfirst!(s::Solver)::Set{Atom}
    top = pop!(s.stack)
    while !isempty(top.free)
        l = Atom(minimum(l->l.I, top.free))
        forked = fork!(top, l)
        if isblocked(top, s.blockers)
            top = isblocked(forked, s.blockers) ? pop!(s.stack) : forked
        elseif !isblocked(forked, s.blockers)
            push!(s.stack, forked)
        end
    end
    top.affirmed
end

function block!(s::Solver, bs::Vector{Blocker})::Nothing
    filter!(s.stack) do a
        !isblocked(a, bs)
    end
    append!(s.blockers, bs)
    nothing
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
    box::Box

    function Onion(rs::Vector{Ribbon}, b::Box)
        @assert length(rs) > 1 "need at least two ribbons, got $length(rs)"
        m = Membrane(GraphDistances(), rs[1])
        @assert isempty(m.targets)
        ds = exterior_distances(m, Set{Atom}())
        new(rs, [Membrane(ds, rs[2])], [Solver(constraints(b, ds))], Set{Atom}[], b)
    end
end

Onion(kernel::Ribbon, basis::Ribbon, n::Int, b::Box = Box(n)) =
    Onion(ribbons(kernel, basis, n), b)

iscomplete(o::Onion) = length(o.diags) == length(o.ribbons)

const counter = Ref(0)

function wrap!(o::Onion, diags::Set{Atom})::Bool
    if length(o.solvers) < length(o.ribbons)
        ds = exterior_distances(o.membranes[end], diags)
        cs = constraints(o.box, ds)
        #println(cs)
        if !issatisfiable(cs)
            #println(violations(cs))
            #println("XX")
            #block!(o.solvers[end], blockers(cs, o.membranes[end].interior, diags))
            return false
        end
        push!(o.solvers, Solver(cs))
        length(o.solvers) < length(o.ribbons) &&
            push!(o.membranes, Membrane(ds, o.ribbons[length(o.solvers) + 1]))
    end
    push!(o.diags, diags)
    true
end

function step!(o::Onion)::Nothing
    iscomplete(o) && pop!(o.diags)

    while true
        while isempty(o.solvers[end])
            length(o.membranes) == length(o.solvers) && pop!(o.membranes)
            pop!(o.solvers)
            isempty(o.diags) && return
            pop!(o.diags)
        end

        counter[] += 1
        diags = popfirst!(o.solvers[end])
        wrap!(o, diags) && break
    end
end

function all_steps!(f, o::Onion)::Nothing
    while !isempty(o.solvers)
        f(o)
        step!(o)
    end
end

function eugrid(o::Onion)
    diags = Matrix{Union{Missing, Bool}}(missing, maximum(o.ribbons[end]).I)
    for (r, d) in zip(o.ribbons, o.diags)
        diags[r] .= in.(r, Ref(d))
    end
    diags
end

function all_eugrids!(o::Onion)::Vector{Matrix{Union{Missing, Bool}}}
    eugrids = Matrix{Union{Missing, Bool}}[]
    all_steps!(o) do o
        push!(eugrids, eugrid(o))
    end
    eugrids
end

function solve!(o::Onion)::Union{BitMatrix, Nothing}
    while !isempty(o.solvers)
        step!(o)
        iscomplete(o) && return eugrid(o)
    end
    nothing
end

solve(args...) = solve!(Onion(args...))

function all_solutions!(o::Onion)::Vector{BitMatrix}
    solutions = BitMatrix[]
    all_steps!(o) do o
        iscomplete(o) && push!(solutions, eugrid(o))
    end
    solutions
end

all_solutions(args...) = all_solutions!(Onion(args...))

function count_solutions(args...)
    n = 0
    i = 0
    all_steps!(Onion(args...)) do o
        iscomplete(o) && (n += 1)
        i += length(o.ribbons)^2
        i % 100_000 == 0 && println(n, " ", state_size.(o.solvers))
    end
    n, i
end

function limit_solve(args...)
    o = Onion(args...)
    counter[] = 0
    k = 1
    while true
        isempty(o.solvers) && return counter[]
        step!(o)
        if iscomplete(o)
            println(counter[])
            return BitMatrix(eugrid(o))
        end
        if counter[] > k
            k *= 2
            println(counter[])
        end
    end
end

#function mem_solve!(o::Onion, ps::Set{BitMatrix})

#6: (4658656, 170719416)
# ~9.9 seconds

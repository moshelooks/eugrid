
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






function Solver(l::Layer)
    isnothing(l.parent) && return Solver([Assignment(Set{Atom}[], l.xx)])



struct Ribbon
end

function wrap!(l::Layer, diags::Set{Atom})
    isnothing(l.child) && return
    for d in diags
        filter!(l.constrained_regions[d]) do r
            delta_distance(l,

        for r in



struct Onion
    borders::Vector{Vector{Atom}}
    layers::Vector{Tuple{Layer, Solver}}

    function Onion(rs::Vector{Ribbon}, b::Box)
        @assert length(rs) > 1 "need at least two ribbons, got $length(rs)"
        ribbiona


        m = Membrane(GraphDistances(), rs[1])
        @assert isempty(m.targets)
        ds = exterior_distances(m, Set{Atom}())
        new(rs, [Membrane(ds, rs[2])], [Solver(constraints(b, ds))], Set{Atom}[], b)
    end
end

function Onion(kernel::Vector{Atom}, basis::Vector{Atom}, n::Int)
    atoms = setdiff!(Set(CartesianIndices((n, n))), kernel)
    borders = [sort(kernel)]
    while !isempty(atoms)
        candidates = map(sum, Iterators.product(borders[end], basis))
        push!(borders, sort!(filter(a->!isnothing(pop!(atoms, a, nothing)), candidates)))
        @assert !isempty(borders[end])
    end
    reverse!(borders)
    Onion(borders, [Layer(n, pop!(borders))])
end




function foo()
    while !isempty(solver)
        diags = solve!(solver)
        isnothing(diags) && return nothing
        constraints = infer_constraints(solver, diags)
        unsatisfiables = unsatisfiable_regions(
        issatisfiable(constraints) && return constraints
        revise!(solver, constraints)
    end
end


            push!(layers, Layer


    for membrane in membranes(kernel, basis, n)
    end
end




function wrap!(o::Onion)
    s = solve!(o.layers[end])
    isnothing(s) && return nothing
    while true
        l = wrap(

function grow!(o::Onion)
    k = length(o.layers)
    while length(o.layers) <= k
        l = wrap!(o)
        if isnothing(l)
            pop!(o.layers)
        else
            push!(o.layers, l)
        end
    end
end

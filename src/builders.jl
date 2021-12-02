struct LayerSolver
    parent::Layer
    constraints::Constraints
    stack::Vector{MonoCNF}
end

function LayerSolver(n::Int, membrane::Vector{Atom})
    parent = Layer(n, membrane)
    ds = membrane_distances(parent, Set{Atom}())
    ts = all_triples(n)
    constraints = Constraints(ds, ts)
    stack = [MonoCNF(constraints)]
    LayerSolver(parent, constraints, stack)
end

function solve!(ls::LayerSolver)
    isempty(ls.stack) && return nothing
    cnf = pop!(stack)

    ds = membrane_distances(ls.parent, diags)
    constraints = Constraints(ds, ls.ts)
    issatisfiable(constraints) && return Layer



                     parent::Layer, diags::Set{Atom})
    ds = membrane_distances(parent, diags)
    constraints = Constraints(ts, ds)







struct OnionBuilder
    membranes::Vector{Vector{Atom}}
    layer_builders::Vector{LayerBuilder}
end

function OnionBuilder(kernel::Vector{Atom}, basis::Vector{Atom}, n::Int)
    atoms = setdiff!(Set(Atoms(n)), kernel)
    membranes = [sort(kernel)]
    while !isempty(atoms)
        candidates = map(sum, Iterators.product(membranes[end], basis))
        push!(membranes, sort!(filter(a->!isnothing(pop!(atoms, a, nothing)), candidates)))
        @assert !isempty(membranes[end])
    end
    reverse!(membranes)
    first_layer = Layer(n, pop!(membranes))
    layer_builders = [LayerBuilder(
    Onion(borders, [Layer(n, pop!(borders))])
end

function grow!(ob::OnionBuilder)
    constraints =

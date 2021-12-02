function constraints(b::Box)
    for r in o.old_regions[n]
        if isfree(parent, r) || freedby(r, diags)
            free!(cs, r)
        elseif freedby(r, diags)
            for d in diags
                if d.a in r && d[r.u] + 1 + foo
                    free!(cs, r)
                    break
                end






struct Layer
    membrane::Membrane
    diags::Set{Atom}
    solver::Solver

    Layer(m::Membrane, diags::Set{Atom}, cs::Constraints) =
        new(m, diags, Solver(
end

Layer(m::Membrane) = Layer(m, Set{Atom}(),
end

function next_layer(l::Layer, m::Membrane)
    while !isempty(l.solver)
        diags = solve!(l.solver)
        ds = exterior_distances(l.membrane, diags)
        constraints = infer_constraints(l.triples, ds)
        issatisfiable(constraints)
        isempty(constraints.violations) && return Layer(m, diags, ds, constraints.clauses)
        revise!(l.solver, constraints.violations)
    end
    nothing
end



regions(t::Triple, a::Atom) = (Region(t, u) for u in max(onexy, a - t.v + onexy):a-onexy)


function wrap!(m::Membrane, diags::Set{Atom})
    ds = exterior_distances(m, diags)
    constraints = infer_constraints(ds)
    unsatisfiable(constraints) && return violations(constraints)
    Membrane(ds

struct LayerSolver
    layer::Layer
    solver::Solver
end

LayerSolver(n::Int, kernel::Vector{Atom}) =
    LayerSolver(Layer(n, kernel), Solver(Set(kernel)))

function wrap!(ls::LayerSolver)::Union{LayerSolver, Nothing}
    diags = solve!(solver)
    isnothing(diags) && return nothing
    ds = membrane_distances(ls.layer, diags)
    cs = Constraints(ls.triples, ds)
    issatisfiable(cs)




Layer(parent::Layer, source_distances::GraphDistances, targets::Vector{Atom}) =
    Layer(merge!(source_distances, parent.frozen), targets)

function target_distances(layer::Layer, diags::Dict{Atom, Bool})::GraphDistances
    ds = GraphDistances()
    for target in layer.targets
        get(layer.frozen, target) do
            tl = target - onexy
            if diags[tl]
                get!(layer.dcache, target) do
                    d = DistanceMatrix(target)
                    ix = CartesianIndices(tl)
                    d.data[ix] .= layer.frozen[tl].data .+ 1
                    d
                end
            else
                t = target - oney
                dt = get(layer.frozen, t) do; ds[t] end
                l = target - onex
                l = get(layer.frozen, l) do; ds[l] end
                get!(layer.ndcache, (objectid(t), objectid(l))) do
                    d = DistanceMatrix(tl)
                    ix = CartesianIndices(tl)
                    d.data[ix] .= min.(view(t.data, ix), view(l.data, ix)) .+ 1
                    d
                end
            end
        end
    end
    ds
end

function region_clauses()
    for (a, ds) in target_distances(layer, diags)
    end
end


function solve!(layer::Layer)::Union{Layer, Nothing}
    while !isempty(layer.cnf)
        diags = popfirst!(layer.cnf)

        for (r,
        ds = target_distances


function eugrid(layer::Layer)
    diags = Matrix{Union{Missing, Bool}}(missing, maximum(keys(layer.frozen)).I)
    for (i, d) in layer.frozen
        minimum(i.I) > 1 && (diags[i - onexy] = isdiag(d))
    end
    diags
end

function extend()



function Layer(parent::Layer, ds::GraphDistances, targets)



    for a in border
        if a in parent.frozen
            x = parent.frozen[a]
        end
        if



foo(layer::Layer, diags::BitVector, ds::GraphDistances, key::DistanceMatrix) = key

foo(layer::Layer, diags::BitVector, ds::GraphDistances, key::Atom) =
    if diags[a]
        get!(layer.dcache, a) do
            r = DistanceMatrix(atom(tl) + onexy)
    br.data[CartesianIndices(tl.data)] .= tl .+ 1
    br
    else
        end
    end
end


function lookup(

function border_distances(layer::Layer, diags::BitVector)::Vector{DistanceMatrix}



    ds = Vector{DistanceMatrix}(undef, length(layer.border))
    for (i, d) in layer.nondeps
        ds[i] = d
    end
    for (i, j, tl) in layer.deps
        if diags[j]
            ds[i] = get!(

function lookup(layer::Layer, solution::Solution{Atom}, l::Int, t::Int, i::Int)
    a, diag = solution[i]
    if diag
        get!(layer.dcache, a)

    diags[i] ? dlookup(

lookup(layer::Layer, ds::DistanceMatrix) = ds
function lookup(layer::Layer, a::Atom)
    if diags[a]
        get(
    else
        nodiag(


    ::DistanceMatrix) = ds


function border_distances(layer::Layer, diags::BitVector)::Vector{DistanceMatrix}
    for ds,

    for (out, tl) in layer.tls
        ds[i] = layer.frozen[tl]
    end
    for (i,

        if diags[tl]
            ds[i] = diag(tl)
        else
            ds[i] = nodiag(tl + onex,

function extend(l::Layer, s::Solution{Atom}, b::Border)::Vector{

struct Extender
    layer::Layer
    solution::

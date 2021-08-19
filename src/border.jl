
struct Distances
    frozen::Dict{Atom, DistanceMatrix}
    dtcache::IdDict{DistanceMatrix, DistanceMatrix}
    dfcache::Dict{NTuple{2, UInt64}, DistanceMatrix}
    alist::LinkedList{Tuple{Atom, DistanceMatrix}}
end

Distances(frozen::Dict{Atom, DistanceMatrix}) =
    Distances(
        frozen,
        IdDict{DistanceMatrix, DistanceMatrix}(),
        Dict{NTuple{2, UInt64}, DistanceMatrix}(),
        nil(Tuple{Atom, DistanceMatrix}))

Distances(n::Int) =
    Distances(Dict(
        i=>DistanceMatrix(reshape(maximum(i.I)-1:-1:0, i.I...))
        for i in Iterators.flatten((CartesianIndices((n, 1)), CartesianIndices((1, n))))))

function freeze(ds::Distances)::Distances
    frozen = copy(ds.frozen)
    for (a, da) in ds.alist
        frozen[a] = da
    end
    Distances(frozen)
end

function dmat(ul::DistanceMatrix)::DistanceMatrix
    n, m = size(ul) .+ 1
    d = DistanceMatrix(undef, n, m)
    d[1:n,m] = n-1:-1:0
    d[n,1:m-1] = m-1:-1:1
    d
end

lookup(ds::Distances, a::Atom)::DistanceMatrix = get(ds.frozen, a) do
    for (i, d) in ds.alist
        i == a && return d
    end
end

function extend(ds::Distances, a::Atom, diag::Bool)::Tuple{Distances, DistanceMatrix}
    if diag
        dul = ds.frozen[a - onexy]
        da = get!(ds.dtcache, dul) do; dmat(dul .+ 1) end
    else
        du = lookup(ds, a - oney)
        dl = lookup(ds, a - onex)
        da = get!(ds.dfcache, (objectid(du), objectid(dl))) do
            ix = CartesianIndices(a.I .- 1)
            dmat(1 .+ min.(view(du, ix), view(dl, ix)))
        end
    end
    Distances(ds.frozen, ds.dtcache, ds.dfcache, cons((a, da), ds.alist))
end

function eugrid(l::Layer)
    diags = Matrix{Union{Missing, Bool}}(missing, maximum(key(l)).I)
    for (br, d) in l.frozen
        minimum(br.I) == 1 && continue
        ul = br - onexy
        diags[ul] = 2 - d[ul]
    end
    for (br, d) in l.alist
        ul = br - onexy
        @assert ismissing(diags[ul])
    = BitMatrix(d[

struct Onion
    eugrid::BitMatrix
    basis::Vector{Atom}
    triples::Vector{Triple}
    layers::Vector{Layer}
end

peel!(o::Onion) = pop!(o.layers)

wrap!(o::Onion) = push!(o.layers, freeze(o.layers[end]))

function eugrid(o::Onion) = BitMatrix(d[



regions(l::Layout, a::Atom) = Iterators.flatten((regions(t, a) for t in l.ts))



function extend(l::Layout, ds::Dict{Atom, DistanceMatrix})::Dict{Atom, DistanceMatrix}
    ds = copy(ds)
end


struct Border
    ds::Dict{Atom, DistanceMatrix}
    active::Set{Atom}
    free::Set{Atom}
    clauses::Vector{Set{Atom}}
end

Border() = Border(Dict(), Set([Atom(1, 1)]), Set(), [])

function extend(layout::Layout, diags::BitMatrix, p::Border=Border())
    ds = copy(p.ds)
    active = Set{Atom}()
    free = Set{Atom}()
    sentinel = Vector{Atom}()
    constraints = Dict{Region, Union{Vector{Atom}, Nothing}}()
    for a in p.active, m in layout.moves
        u = a + m
        haskey(ds, u) && continue
        push!(active, u)
        ds[u] = du = DistanceMatrix(undef, u.I)
        v = u - onexy
        du[:, u[2]] .= v[1]:-1:0
        du[u[1], :] .= v[2]:-1:0
        minimum(u.I) == 1 && continue
        ix = CartesianIndices(v.I)
        @views du[ix] = 1 .+ (diags[v] ? ds[v][ix] : min.(ds[v+onex][ix], ds[v+oney][ix]))

        negated = false
        clausal_regions = Vector{Region}()
        for r in regions(layout, u)
            if negates(r, du)
                negated = true
                constraints[r] = nothing
            elseif !isnothing(get(constraints, sentinel))
                bound = affirmation_bound(r, du)
                if bound < 0
                    constraints[r] = nothing
                elseif bound == 0 && !negated
                    push!(clausal_regions, r)
                elseif !haskey(constraints, r)
                    constraints[r] = Vector{Int}()
                end
            end
        end
        if negated
            for r in clausal_regions
                if !haskey(constraints, r)
                    constraints[r] = Vector{Int}()
                end
            end
        else
            push!(free, u)
            for r in clausal_regions
                if haskey(constraints, r)
                    push!(constraints[r], u)
                else
                    constraints[r] = [u]
                end
            end
        end
    end

    clauses = Vector{Set{Atom}}()
    for (r, c) in constraints
        isnothing(c) && continue
        isempty(c) && return nothing
        push!(clauses, Set(c))
    end
    Border(ds, active, free, clauses)
end

function grow()

#=

        for r in regions(u)
            negated |= update_clause!(


    end
    Border(ds, active)
end

function atoms(b::Border, r::Region) =

function constraints(ts::Vector{Triple}, b::Border)


function pose(ts::Vector{Triple}, b::Border)::Union{MonoCNF, Nothing}
    clauses = Vector{Vector{Atom}}()
    free = Set{Atom}()
    for (a, constraints) in atom_constraints
        negated = false
        for c in constraints
            if update_clause!(c, b.ds[a], nothing)
                negated = true
                p
        end
        !negated && push!(a, free)
    end
    for (c, atoms) in constraint_atoms
        for a in atoms
            update_clause!(c, b.ds[a], !(a in free))


    constraints = Dict{Tuple{Triple, Atom}, Vector{Atom
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
=#

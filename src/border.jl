const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

struct Triple
    a::Int
    b::Int
    c::Int
    Triple(a::Int, b::Int) = new(a, b, sqrt(a^2 + b^2))
end

const Atom = CartesianIndex{2}

struct Region
    t::Triple
    u::Atom
end

regions(t::Triple, v::Atom) =
    (Region(t, u) for u in max(onexy, v - CartesianIndex(t.a + 1, t.b + 1)):v-onexy)

struct Layout
    ts::Vector{Triple}
    moves::Vector{Atom}
end

regions(l::Layout, v::Atom) =
    Iterators.flatten((regions(t, v) for t in l.ts))

const DistanceMatrix = Matrix{Int}

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

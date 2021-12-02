


struct Solution
    clauses::Vector{Set{Atom}}
    free::Set{Atom}
    diags::Set{Atom}
end

Solution(clauses, free) =
    simplify!(Solution([Set{Atom}(c) for c in clauses], Set{Atom}(free), Set{Atom}()))

iscomplete(solution::Solution)::Bool = isempty(solution.free)

function simplify!(solution::Solution)::Solution
    unique!(solution.clauses)

    for subset in sort(solution.clauses, by=length)
        filter!(solution.clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end

    filter!(solution.clauses) do c
        length(c) > 1 && return true
        l = only(c)
        pop!(solution.free, l)
        push!(solution.diags, l)
        false
    end

    solution
end

function fork!(score, solution::Solution)::NTuple{2, Solution}
    _,



struct Solver
    layer::Layer
    clauses::Dict{Region, Vector{Atom}}
    stack::Vector{Solution}

    Solver(layer, clauses) = new(layer, clauses, [Solution(clauses)])
end



const Violations = Vector{Region}

function Base.isempty(solver::Solver)::Bool = isempty(solver.stack)

struct ClauseBuilder
    layout::Layout

end

const

function constraints(o::Onion, ds::GraphDistances)::Dict{Region, Union{Vector{Atom}, Nothing}}

function extend(ds::GraphDistances, ts::Vector{Triple}, diags::Vector{Tuple{Atom, Bool}})
    e = Extender(ds)
    b = Builder()
    for (a, diag) in diags
        e, da = extend(e, a, diag)
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


    (ds::GraphDistances

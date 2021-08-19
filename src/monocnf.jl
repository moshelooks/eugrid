const Clauses{T} = Vector{Set{T}}

function simplify!(clauses::Clauses{T})::Clauses{T}
    unique!(clauses)
    for subset in sort(clauses, by=length)
        filter!(clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end
    clauses
end

const Assignment{T} = Dict{T, Bool}

struct Solution{T}
    clauses::Clauses{T}
    assignment::Assignment{T}
    free::Set{T}
end

function fork!(s::Solution{T}, l::T)::NTuple{2, Solution{T}}
    pop!(s.free, l)

    affirmed = Solution(
        [copy(c) for c in s.clauses if !(l in c)],
        push!(copy(s.assignment), l=>true),
        copy(s.free))

    for c in s.clauses
        pop!(c, l, nothing)
    end
    affirm_singletons!(s)

    affirmed, s
end



struct MonoCNF{T}
    domain::Vector{T}
    clauses::Vector{Set{T}}
    stack::Vector{Solution{T}}
    affirmed::Vector{T}

    function MonoCNF(clauses=Vector{Set{Int}}(), free=Set{Int}(), affirmed=Vector{Int}())
        cnf = new(clauses, free, affirmed)
        validate(cnf)
        return cnf
    end
end

function solve!(cnf::MonoCNF{T})::Union{Solution{T}, Nothing}
    n = length(cnf.cnf)
    while !isempty(cnf.stack)
        top = pop!(cnf.stack)
        length(top) == n && return top


function validate(cnf::MonoCNF)::Nothing
    for c in cnf.clauses
        for l in c
            @assert l in cnf.free "free missing $l in $cnf"
        end
    end
    for a in cnf.affirmed
        @assert !(a in cnf.free) "$a in free and affirmed in $cnf"
    end
end

function solve!(cnf::MonoCNF{T}, heuristic=nothing)::Union{Vector{Tuple{T, Bool}}, Nothing}
    isempty(cnf.stack) && return nothing
    if !isnothing(heuristic)
        scores = Iterators.map(stack) do
        i = argma


struct Solver
end

function


function affirm_singletons!(cnf::MonoCNF)::MonoCNF
    simplify!(cnf)
    filter!(cnf.clauses) do c
        length(c) > 1 && return true
        l = only(c)
        pop!(cnf.free, l)
        push!(cnf.affirmed, l)
        false
    end
    cnf
end

function simplify!(cnf::MonoCNF)::MonoCNF
    unique!(cnf.clauses)
    for subset in sort(cnf.clauses, by=length)
        filter!(cnf.clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end
    cnf
end

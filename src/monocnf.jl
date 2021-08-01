struct MonoCNF
    clauses::Vector{Set{Int}}
    free::Set{Int}
    affirmed::Vector{Int}

    function MonoCNF(clauses=Vector{Set{Int}}(), free=Set{Int}(), affirmed=Vector{Int}())
        cnf = new(clauses, free, affirmed)
        validate(cnf)
        return cnf
    end
end

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

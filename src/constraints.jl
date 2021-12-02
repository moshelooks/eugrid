struct Constraints
    domain::Set{Atom}
    regions::Dict{Region, Union{Vector{Atom}, Nothing}}
end

const _sentinel = [Atom(0,0)]

free!(cs::Constraints, r::Region) = setindex!(cs.regions, r, nothing)
isfree(cs::Constraints, r::Region)::Bool = isnothing(get(cs.regions, r, _sentinel))

cover!(c::Constraints, r::Region) = get!(Vector{Atom}, c.regions, r)

function infer!(constraints::Constraints, ts::Vector{Triple}, d::DistanceMatrix)::Nothing
    negated = false
    clausal_regions = Vector{Region}()
    for r in regions(ts, d.a)
        delta = distance_delta(d, r)
        if delta == negation_delta(r, u)
            negated = true
            free!(constraints, r)
        elseif !isfree(contraints, r)
            delta -= affirmation_delta(r, u)
            if delta < 0
                free!(constraints, r)
            elseif delta == 0 && !negated
                push!(clausal_regions, r)
            else
                cover!(c, r)
            end
        end
    end
    if negated
        foreach(clausal_regions) do r; cover!(c, r) end
    else
        push!(c.domain, d.a)
        foreach(clausal_regions) do r; push!(cover!(c, r), d.a) end
    end
end

function MonoCNF(clauses::Vector{Set{Int}})::MonoCNF
    domain = unique(l for c in clauses for l in c)
    counts = Dict(l=>0 for l in domain)
    for c in clauses, l in c
        counts[l] += 1
    end
    MonoCNF(clauses, counts)
end

function true!(cnf::MonoCNF, l::Int)

singletons(cnf:: MonoCNF)::Vector{Int} = unique(only(c) for c in cnf if length(c) == 1)


function extract_singletons!(cnf::MonoCNF)::Vector{Int}
    extracted = singletons(cnf)
    filter!(c->!any(in(c), toremove), cnf)
    end
end

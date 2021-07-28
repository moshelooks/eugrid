module Pythagorean

export Triple, primitive_triples, scaled_triples

const ABC = [1 -2 2; 2 -1 2; 2 -2 3], [1 2 2; 2 1 2; 2 2 3], [-1 2 2; -2 1 2; -2 2 3]

children(x::Vector{Int}) = [m * x for m in ABC]

Triple = Tuple{Int, Int, Int}

function primitive_triples(n::Int)::Vector{Triple}
    triples = Vector{Triple}()
    stack = [[3, 4, 5]]
    while !isempty(stack)
        triple = pop!(stack)
        a, b = minmax(triple[1:2]...)
        a + b > n && continue
        push!(triples, (a, b, triple[3]))
        push!(stack, children(triple)...)
    end
    sort!(triples, by=x->(x[1]+x[2]))
end

function scaled_triples(n::Int)::Vector{Triple}
    triples = Vector{Triple}()
    for triple in primitive_triples(n)
        a, b, c = triple .* div(n, sum(triple[1:2]))
        push!(triples, (a, b, c), (b, a, c))
    end
    triples
end

function span(t::Triple)
    a, b = minmax(t[1], t[2])
    a + div(b - a, 2)
end

end # module

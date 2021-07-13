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
        a > n && continue
        b <= n && push!(triples, (a, b, triple[3]))
        push!(stack, children(triple)...)
    end
    sort!(triples, by=x->(x[3], x[1]))
end

function scaled_triples(n::Int)::Vector{Triple}
    triples = Vector{Triple}()
    for triple in primitive_triples(n)
        a, b, c = triple .* div(n, triple[2])
        push!(triples, (a, b, c), (b, a, c))
    end
    triples
end

end # module

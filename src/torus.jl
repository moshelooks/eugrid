struct Torus
    n::Int
    k::Int
    d::Array{Int, 4}
    diags::CircularArray{Bool, 2, BitMatrix}
    triples::Matrix{Vector{Tuple{CartesianIndex{2}, Int}}}
end

function Torus(n::Int, k::Int=n)::Torus
    d = repeat(CartesianIndices((k, k)) .|> x->sum(x.I), outer=(1, 1, n, n))
    triples = [Vector{Tuple{CartesianIndex{2}, Int}}() for _ in onexy:(onexy * k)]
    for (a, b, c) in Pythagorean.scaled_triples(k)
        triple = CartesianIndex(a, b), c
        for j in onexy:triple[1]
            push!(triples[j], triple)
        end
    end
    Torus(n, k, d, CircularArray(falses(n, n)), triples)
end

Base.CartesianIndices(t::Torus) = CartesianIndices(t.diags)

wrap(t::Torus, i::CartesianIndex{2})::CartesianIndex =
    CartesianIndex(mod.(i.I, axes(t.diags)))

deltas(t::Torus)::Array{Int, 3} =
    [t.d[u, i] - l for i in CartesianIndices(t), (u, l) in t.triples[onexy]]

bound(t::Torus)::Int = maximum(abs.(deltas(t)))

dvia(t::Torus, u::CartesianIndex{2}, v::CartesianIndex{2}, w::CartesianIndex{2})::Int =
    (minimum(v.I) == 1 ? maximum(v.I) : t.d[v - onexy, u] + 1) +
    (minimum(w.I) == 0 ? maximum(w.I) : t.d[w, wrap(t, u + v)])

function add_diag(t::Torus, u::CartesianIndex{2})::Nothing
    @assert !t.diags[u]
    t.diags[u] = true

    du = view(t.d, :, :, u)
    du[:, 1] = du[1, :] = 1:t.k
    diags = view(t.diags, u+onexy:u+(t.k - 1)*onexy)
    for i in CartesianIndices(diags)
        du[i+onexy] = 1 + (diags[i] ? du[i] : min(du[i+onex], du[i+oney]))
    end

    kxy = t.k * onexy
    for i in Iterators.take(CartesianIndices(du), t.k^2 - 1)
        v = wrap(t, u - kxy + i)
        dvu = maximum(i.I) == t.k ? t.k - minimum(i.I) : t.d[kxy - i, v]
        dv = view(t.d, onexy + kxy - i: kxy, v)
        dv[1] == dvu + 1 && continue
        dv .= min.(dv, dvu .+ view(du, onexy:i))
    end
    nothing
end

symmetrize(t::Torus, u::CartesianIndex{2})::BitMatrix =
    t.diags[u:u+onexy*(t.k-1)].data
    #BitMatrix(t.diags[i + u - onexy] for i in CartesianIndices(t))

#=
function symmetrize(t::Torus, u::CartesianIndex{2})::BitMatrix
    sym = BitMatrix(t.diags[i + u - onexy] for i in CartesianIndices(t))
    sym2 = transpose(sym)
    for (a, b) in zip(sym, sym2)
        if a
            !b && (return sym2)
        elseif b
            break
        end
    end
    sym
end
=#

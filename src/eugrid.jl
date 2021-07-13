module eugrid

using CircularArrays: CircularArray

const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

TripleVec = Vector{Tuple{CartesianIndex{2}, Int}}

function pythagorean_triples(n)::TripleVec
    ts = TripleVec()
    for (a, b, c) in [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25),
                      (20, 21, 29), (12, 35, 37), (9, 40, 41), (28, 45, 53),
                      (11, 60, 61), (16, 63, 65), (33, 56, 65), (48, 55, 73),
                      (13, 84, 85), (36, 77, 85), (39, 80, 89), (65, 72, 97)]
        b > n && continue
        a, b, c = [a, b, c] .* div(n, b)
        push!(ts, (CartesianIndex(a, b), c), (CartesianIndex(b, a), c))
    end
    ts
end

struct Torus
    n::Int
    k::Int
    d::Array{Int, 4}
    diags::CircularArray{Bool, 2, BitMatrix}
    triples::Matrix{TripleVec}
end

Base.CartesianIndices(t::Torus) = CartesianIndices(t.diags)

function Torus(n::Int, k::Int=n)::Torus
    d = repeat(CartesianIndices((k, k)) .|> x->sum(x.I), outer=(1, 1, n, n))
    triples = [TripleVec() for _ in onexy:(onexy * k)]
    for (u, l) in pythagorean_triples(n)
        for j in onexy:u
            push!(triples[j], (u, l))
        end
    end
    Torus(n, k, d, CircularArray(falses(n, n)), triples)
end

wrap(t::Torus, i::CartesianIndex{2})::CartesianIndex =
    CartesianIndex(mod.(i.I, axes(t.diags)))

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

deltas(t::Torus) = [t.d[u, i] - l for i in CartesianIndices(t), (u, l) in t.triples[onexy]]

bound(t::Torus)::Int = maximum(abs.(deltas(t)))

function dvia(t::Torus,
              u::CartesianIndex{2}, v::CartesianIndex{2}, w::CartesianIndex{2})::Int
    duv = minimum(v.I) == 1 ? maximum(v.I) : t.d[v - onexy, u] + 1
    dvw = minimum(w.I) == 0 ? maximum(w.I) : t.d[w, wrap(t, u + v)]
    duv + dvw
end

function addable(t::Torus, b::Int, u::CartesianIndex{2})::Bool
    t.diags[u] && return false
    for j in onexy:(onexy * t.k)
        i = wrap(t, u - j + onexy)
        for (k, l) in t.triples[j]
            d_i_k = t.d[k, i] - 1
            d_i_k >= l - b  && continue
            dvia(t, i, j, k - j) == d_i_k && return false
        end
    end
    true
end

function score(t::Torus, b::Int)::Matrix{Int}
    scores = zeros(Int, t.n, t.n)
    for u in CartesianIndices(t)
        !addable(t,b,u) && continue
        c = deepcopy(t)
        add_diag(c, u)
        scores[u] = count(addable(c, b, i) for i in CartesianIndices(t))
    end
    scores
end

function tdist(n::Int, u::CartesianIndex{2}, v::CartesianIndex{2})::Int
    x = abs.(v.I .- u.I)
    maximum(min.(x, n .- x))
end

function add_all2(t::Torus, b)::Nothing
    added = []
    while true
        scores = score(t, b)
        s = maximum(scores)
        s < 1 && return
        cs = [i for i in CartesianIndices(t) if scores[i] == s]
        for a in reverse(added)
            length(cs) == 1 && break
            tds = [tdist(t.n, c, a) for c in cs]
            mtd = maximum(tds)
            cs = [c for (c, td) in zip(cs,tds) if td == mtd]
        end
        ix = cs[1]
        println(s, " ", ix)
        add_diag(t, ix)
        push!(added, ix)
    end
    nothing
end

function add_all(t::Torus)::Nothing
    for b in -1:1
        println("xxx ", b)
        add_all2(t, b)
    end
    nothing
end

function exnil(n::Int, k::Int=n)::Torus
    t = Torus(n, k)
    for i in CartesianIndices(t)
        sum(i.I.%2) == 0 && add_diag(t, i)
    end
    add_all(t)
    t
end

end # module

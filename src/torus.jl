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
    @assert !t.diags[u] u
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

using Evolutionary
using Statistics

function frombits(xs)
    k = Int(sqrt(length(xs)))
    t = Torus(k, k)
    ix = CartesianIndices(t.diags)
    for i in eachindex(xs)
        if xs[i]
            add_diag(t, ix[i])
        end
    end
    t
end


#maxabs(t) = maximum(abs.(err(t)))
#meanabs(t) = mean(abs.(err(t)))
#mse(t) = mean(err(t).^2)
#errstats(t) = (maxabs(t), meanabs(t), mse(t))

function oscore(xs)
    t = frombits(xs)
    #se = ((t.d[u, i] - l)^2 for i in CartesianIndices(t), (u, l) in t.triples[onexy])
    se = 0
    for j in 1:t.k, k in 1:t.k
        de = sqrt(j^2+k^2)
        for i in CartesianIndices(t)
            se += (t.d[j, k, i] - de)^2
        end
    end

    1e8 * bound(t) + mean(se) #* 1e2 + sum(xs)
    #div(1e6 * maxabs(t), 1) * 1000 + meanabs(t)
    #mse(t)
end

function ga(k, inst=falses(k^2))
    Evolutionary.optimize(oscore, inst, GA(
        populationSize=40000,
        crossoverRate=0.5,
        mutationRate=0.1,
        epsilon=0.1,
        selection=susinv,
        crossover=discrete,
        mutation=flip), Evolutionary.Options(iterations=100, show_trace=true))
end



function frombits2(xs)
    n = Int(length(xs)/2)
    t = Torus(n, n)
    for i in 1:n
        !xs[i] && continue
        for j in 1:n
            u = wrap(t, CartesianIndex(n-j+1,i+j-1))
            add_diag(t, u)
        end
    end
    #=
    for i in (n+1):(2*n)
        !xs[i] && continue
        for j in 1:n
            u = CartesianIndex(j, i-n)
            t.diags[u] && continue
            #t.diags[wrap(t, u - onexy)] && continue
            add_diag(t, u)
        end
    end
    for i in (2*n+1):(3*n)
        !xs[i] && continue
        for j in 1:n
            u = CartesianIndex(i - 2*n, j)
            t.diags[u] && continue
            #t.diags[wrap(t, u - onexy)] && continue
            add_diag(t, u)
        end
    end
    =#
    for i in (n+1):(2*n)
        !xs[i] && continue
        k = i - n
        for x in k:k:n
            for y in k:k:n
                u = CartesianIndex(x, y)
                t.diags[u] && continue
                add_diag(t, u)
            end
        end
    end
    t
end


function oscore2(xs)
    t = frombits2(xs)
    se = ((t.d[u, i] - l)^2 for i in CartesianIndices(t), (u, l) in t.triples[onexy])
    #=
    se = 0
    for j in 1:t.k, k in 1:t.k
        de = sqrt(j^2+k^2)
        for i in CartesianIndices(t)
            se += (t.d[j, k, i] - de)^2
        end
    end
    =#
    1e8 * bound(t) + mean(se) #* 1e2 + sum(xs)
    #div(1e6 * maxabs(t), 1) * 1000 + meanabs(t)
    #mse(t)
end

function ga2(k, inst=falses(2*k))
    Evolutionary.optimize(oscore2, inst, GA(
        populationSize=1000,
        crossoverRate=0.5,
        mutationRate=0.05,
        epsilon=0.1,
        selection=susinv,
        crossover=discrete,
        mutation=flip), Evolutionary.Options(iterations=100, show_trace=true))
end

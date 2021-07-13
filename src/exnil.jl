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

struct Score
    t::Torus
    u::CartesianIndex{2}
    addables::Vector{CartesianIndex{2}}
end

Score(t::Torus, b::Int)::Score =
    Score(t, CartesianIndex(0, 0), filter(u->addable(t,b,u), CartesianIndices(t)))

function score(prev::Score, b::Int, u::CartesianIndex, best::Int)
    t = deepcopy(prev.t)
    add_diag(t, u)
    addables = Vector{CartesianIndex{2}}()
    for (n, i) in enumerate(prev.addables)
        addable(t, b, i) && push!(addables, i)
        length(addables) + length(prev.addables) - n <= best + 1 && return nothing, best
    end
    Score(t, u, addables), length(addables)
end


function add_all2(t::Torus, b::Int)::Torus
    state = Score(t, b)
    while true
        next, best = state, 0
        symmetries = Set{BitMatrix}()
        for u in state.addables
            sym = symmetrize(state.t, u)
            sym in symmetries && continue
            transpose(sym) in symmetries && continue
            push!(symmetries, sym)
            next_u, best = score(state, b, u, best)
            if !isnothing(next_u)
                next = next_u
                best == length(state.addables) -1 && break
            end
        end
        best < 1 && break
        state = next
        println(best, " ", state.u)
    end
    state.t
end

function add_all(t::Torus)::Torus
    for b in -1:1
        println("xxx ", b)
        t = add_all2(t, b)
    end
    t
end

function exnil(n::Int, k::Int=n)::Torus
    t = Torus(n, k)
    for i in CartesianIndices(t)
        sum(i.I.%2) == 0 && add_diag(t, i)
    end
    add_all(t)
end

#=

struct Adder
    t::Torus
    b::Int
    addables::Vector{CartesianIndex{2}}
end

Adder(t::Torus, b::Int)::Adder = Adder(t, b, filter(u->addable(t,b,u), CartesianIndices(t)))

function add_diag(t::Torus, u::CartesianIndex{2}, a::Adder)::Nothing
    add_diag(t, u)
    filter!(u->addable(t,b,u), a.addables)
    nothing
end

function score(a::Addable)::Matrix{Int}
    scores = zeros(Int, a.t.n, a.t.n)

    kxy = t.k * onexy
    for u in a.addables
        du = t.d[ :, :, u]
        du[:, 1] = du[1, :] = 1:t.k
        diags = view(t.diags, u+onexy:u+(t.k - 1)*onexy)
        for i in CartesianIndices(diags)
            du[i+onexy] = 1 + (diags[i] ? du[i] : min(du[i+onex], du[i+oney]))
        end

        dvs = [du]

        for i in Iterators.take(CartesianIndices(du), t.k^2 - 1)
            v = wrap(t, u - kxy + i)
            dvu = maximum(i.I) == t.k ? t.k - minimum(i.I) : t.d[kxy - i, v]
            dv = t.d[onexy + kxy - i: kxy, v]
            dv[1] == dvu + 1 && continue
            dv .= min.(dv, dvu .+ view(du, onexy:i))
            push!(dvs, dv)
        end

        for i in CartesianIndices(du)
            v = wrap(t, u - kxy + i)



        c = deepcopy(t)
        add_diag(c, u)
        scores[u] = count(addable(c, b, i) for i in a.addables)
    end
    scores
end
=#

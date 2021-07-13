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

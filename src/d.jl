function d(t::Torus, u::CartesianIndex{2}, v::CartesianIndex{2})::Int
    duv = 0
    while true
        v[1] == 0 && return duv+v[2]
        v[2] == 0 && return duv+v[1]
        !t.diags[u] && break
        u = wrap(t, u + onexy)
        v -= onexy
        duv += 1
    end
    while t.diags[u+v-onexy]
        v -= onexy
        duv += 1
        v[1] == 0 && return duv+v[2]
        v[2] == 0 && return duv+v[1]
    end
    while true
        duv += 1
        if v[1] > t.rbound[v[2], u]
            u += onex
            v -= onex
            v[1] == 0 && return duv+v[2]
        else
            u += oney
            v -= oney
            v[2] == 0 && return duv+v[1]
        end
    end
end


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

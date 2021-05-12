module eugrid

export shortest_paths, Grid

import LinearAlgebra as la

function shortest_paths(diags::BitMatrix)::Tuple{Matrix{Int},Matrix{Set{CartesianIndex{2}}}}
    n, m = size(diags)
    d = (0:n) .+ (0:m)'
    b = similar(d, Set{CartesianIndex{2}})
    b[:, 1] .= (eltype(b)(),)
    b[1, :] .= (b[1, 1],)

    for i = 1:n
        for j = 1:m
            shortest = min(d[i+1, j], d[i, j+1])
            if diags[i, j] && d[i, j] <= shortest
                if d[i, j] == shortest
                    b[i+1, j+1] = intersect(b[i+1, j], b[i, j], b[i, j+1])
                elseif d[i, j] <= shortest
                    shortest = d[i, j]
                    b[i+1, j+1] = union(b[i, j], (CartesianIndex(i, j),))
                end
            elseif d[i+1, j] == d[i, j+1]
                b[i+1, j+1] = intersect(b[i+1, j], b[i, j+1])
            elseif d[i+1, j] == shortest
                b[i+1, j+1] = b[i+1, j]
            else
                b[i+1, j+1] = b[i, j+1]
            end
            d[i+1, j+1] = shortest + 1
        end
    end
    d, b
end


struct Grid
    n_inner::Int
    diags::BitMatrix

    dil::Array{Int16,3}
    dit::Array{Int16,3}

    dg::Matrix{Int}
    de::Matrix{Float64}

    deltas::Matrix{Float64}
    scores::Matrix{Float64}

    function Grid(n::Int)
        n_inner = n - 1
        diags = falses(n_inner, n_inner)

        iix() = CartesianIndices((1:n_inner, 0:n_inner-1)) .|> x -> x.I

        #dil = [la.norm(inner .- outer, 1) for inner in iix(), outer in lindices(n_inner)]
        dil = [convert(Int16, la.norm(inner .- outer, 1)) for inner in iix(), outer in lindices(n_inner)]
        dit = permutedims(dil, (2, 1, 3))

        dg = [la.norm(l .- t, 1) for l in lindices(n_inner), t in tindices(n_inner)]
        de = [la.norm(l .- t, 2) for l in lindices(n_inner), t in tindices(n_inner)]

        deltas = -2 .* (de .- dg) .- 1
        scores = Array{Float64,2}(undef, n_inner, n_inner)
        #=
        Threads.@threads for ij in  CartesianIndices(scores)
            i, j = ij.I
            #delta = view(deltas, i:n_inner+j-1, j:n_inner+i-1)
            #scores[ij] = sum(delta)
            scores[ij] = sum(
                deltas[d] for d in CartesianIndices((i:n_inner+j-1, j:n_inner+i-1)))
        end
        =#

        new(n_inner, diags, dil, dit, dg, de, deltas, scores)
    end

end

lindices(n_inner::Int) =
    zip([1:n_inner; fill(n_inner, n_inner - 1)], [fill(0, n_inner); 1:n_inner-1])

tindices(n_inner::Int) = reverse.(lindices(n_inner))

function add(g::Grid, x::Int, y::Int)::Nothing
    @assert !g.diags[x, y]
    g.diags[x, y] = true

    lix = x:g.n_inner+y-1
    dxyl = g.dil[x, y, lix]

    if x > 1 && y < g.n_inner
        sps, bb = shortest_paths(g.diags[x-1:-1:2, y+1:g.n_inner-1])
        dil = view(g.dil, x-1:-1:1, y+1:g.n_inner, lix)
        dil .= min.(dil, reshape(dxyl, 1, 1, :) .+ sps .+ 1)
    end

    tix = y:g.n_inner+x-1
    dxyt = g.dit[x, y, tix]

    if x < g.n_inner && y > 1
        sps, bb = shortest_paths(g.diags[x+1:g.n_inner-1, y-1:-1:2])
        dit = view(g.dit, x+1:g.n_inner, y-1:-1:1, tix)
        dit .= min.(dit, reshape(dxyt, 1, 1, :) .+ sps)
    end

    dg = view(g.dg, lix, tix)
    dg .= min.(dg, reshape(dxyl, :, 1) .+ reshape(dxyt, 1, :) .+ 1)

    return
end

end # module

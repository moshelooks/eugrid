abstract type GrowthHeuristic end

function line_score(delta2::Int, tl::Int, br::Int)::Float64
    (tl += 1) == br && return 0.0
    delta = sqrt(delta2)
    abs(delta - br) - abs(delta - tl)
end

function grow_corner_diags(n::Int, leq=false)::BitMatrix
    ds = Matrix{Int}(undef, n+1, n+1)
    ds[:, 1] = ds[1, :] = 0:n
    diags = BitMatrix(undef, n, n)
    pred = leq ? >=(0) : >(0)
    for i in CartesianIndices(diags)
        tl = ds[i]
        l = ds[i+onex]
        t = ds[i+oney]
        br = min(l, t) + 1
        if pred(line_score(sum(i.I.^2), tl, br))
            ds[i + onexy] = tl + 1
            diags[i] = true
        else
            ds[i + onexy] = br
            diags[i] = false
        end
    end
    diags
end

gamma_weight(opp, adj, x::Int, y::Int, ::Val{:unweighted})::Float64 = 1.0
gamma_weight(opp, adj, x::Int, y::Int, ::Val{:angular})::Float64 = atan(opp / adj)
gamma_weight(opp, adj, x::Int, y::Int, ::Val{:minmax})::Float64 = /(minmax(x, y)...)^-1
gamma_weight(opp, adj, x::Int, y::Int, ::Nothing=nothing)::Float64 =
    gamma_weight(opp, adj, x, y, Val(:angular)) * gamma_weight(opp, adj, x, y, Val(:minmax))

function gamma_score(v::Vertex, tl, br, W=nothing)::Float64
    vsum = sum(v.I)
    v2 = v.I.^2
    total_score = 0.0

    for i in 1:v[1]-1
        delta2 = i^2 + v2[2]
        s = line_score(delta2, tl[i], br[i])
        s != 0 && (total_score += gamma_weight(v[2], delta2 - 0.25, i, v[2], W) * s)
    end

    delta2 = sum(v2)
    s = line_score(delta2, tl[v[1]], br[v[1]])
    s != 0 && (total_score += gamma_weight(vsum - 0.5, 2 * delta2 - vsum, v.I..., W) * s)

    for i in v[1]+1:vsum-1
        j = vsum - i
        delta2 = v2[1] + j^2
        s = line_score(delta2, tl[i], br[i])
        s != 0 && (total_score += gamma_weight(v[1], delta2 - 0.25, v[1], j, W) * s)
    end

    total_score
end

function propagate_distances!(i::Int, parents, children)
    n = size(parents, 1)
    l = view(parents, :, i)
    t = view(parents, :, i+1)
    br = view(children, :, i)
    br[1] = l[1] + 1
    br[n+1] = t[n] + 1
    view(br, 2:n) .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1
end

#function diag!(tl, br)
#    br .= view(tl, 1


function gamma_score!(v::Vertex, grandparents, parents, children, avoid=nothing; W=nothing)
    scores = Vector{Float64}(undef, size(children, 2))
    Threads.@threads for i in 1:size(children, 2)
        v_i = v + CartesianIndex(-1, 1) * (i - 1)
        tl = view(grandparents, :, i)
        br = propagate_distances!(i, parents, children)
        scores[i] = avoid[v_i] ? -1 : gamma_score(v_i, tl, br, W)
    end
    scores
end

function sparsity_cutoff(scores, sparsity)
    sparsity >= 1 && return maximum(scores) + 1
    partialsort(scores, 1 + floor(Int, sparsity * length(scores)))
end

function generate!(diags, i::Int, grandparents, parents, children, avoid=nothing;
                   W=nothing, sparsity=nothing)
    n = checksquare(diags)
    depth = size(parents, 1) - 1
    if i <= n
        v_0 = Vertex(i, 1)
        width = size(children, 2) - 2
        children[1:depth+2, 1] .= 0:depth+1
        children[1:depth+2, width+2] .= depth+1:-1:0
        children = view(children, :, 2:width+1)
    else
        v_0 = Vertex(n, i-n+1)
        width = size(children, 2)
    end
    scores = gamma_score!(v_0, grandparents, parents, children, avoid, W=W)
    cutoff = isnothing(sparsity) ? 0.0 : max(0.0, sparsity_cutoff(scores, sparsity))
    @Threads.threads for i in 1:width
        v_i = v_0 + CartesianIndex(-1, 1) * (i - 1)
        if scores[i] >= cutoff
            diags[v_i] = true
            children[2:depth+1, i] .= view(grandparents, 1:depth, i) .+ 1
        else
            diags[v_i] = false
        end
    end
end

const Buffer = Array{Int, 3}

diags_buffer(n::Int) = Buffer(undef, 2*n+1, n+2, 3)

order(buffer::Buffer)::Int = size(buffer, 2) -2

function children_of(buffer::Buffer, i::Int)
    n = order(buffer)
    view(buffer, 1:i+2, 1:(i <= n ? i+2 : 2n-i), mod1(i+2, 3))
end

function next_grandparents_of(buffer::Buffer, i::Int, parents)
    n = order(buffer)
    i >= n ? view(parents, :, 2:(i == n ? n : 2n+1-i)) : parents
end

function next_parents_of(buffer::Buffer, i::Int, children)
    n = order(buffer)
    i == n ? view(children, :, 2:n+1) : children
end

function grow_diags!(buffer::Buffer, diags, avoid; W=nothing, sparsity=nothing)
    n = order(buffer)

    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:2n-1
        children = children_of(buffer, i)
        generate!(diags, i, grandparents, parents, children, avoid, W=W, sparsity=sparsity)
        grandparents = next_grandparents_of(buffer, i, parents)
        parents = next_parents_of(buffer, i, children)
    end

    diags
end

function grow_grid(n::Int; W=nothing, sparsity=nothing, margin=n)
    m = n + margin
    avoid = falses(m, m)
    buffer = diags_buffer(m)
    base = grow_diags!(buffer, BitMatrix(undef, m, m), avoid, W=W, sparsity=sparsity)
    diags = base[(margin+1)*onexy:m*onexy]
    #=
    avoid = falses(m, m)
    avoid[(margin+1)*onexy:m*onexy] .= view(diags, n:-1:1, :)
    grow_lower!(buffer, base, avoid, W=W, sparsity=sparsity))
    antidiags = base[(n+1)*onexy:2n*onexy]

    fixme grow lower assumes lower generalize grow_Base!? masked array?
arctang.ent=#
end

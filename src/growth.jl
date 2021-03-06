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
gamma_weight(opp, adj, x::Int, y::Int, ::Val{:angular})::Float64 = atan(opp, adj)
gamma_weight(opp, adj, x::Int, y::Int, ::Val{:minmax})::Float64 = /(minmax(x, y)...)^-1
gamma_weight(opp, adj, x::Int, y::Int, ::Nothing=nothing)::Float64 =
    gamma_weight(opp, adj, x, y, Val(:angular)) * gamma_weight(opp, adj, x, y, Val(:minmax))

function gamma_score(v::Vertex, tl, br, W=nothing)::Float64
    vsum = sum(v.I)
    v2 = v.I.^2
    summands = Float64[]

    for i in 1:v[1]-1
        delta2 = i^2 + v2[2]
        s = line_score(delta2, tl[i], br[i])
        s != 0 && push!(summands, gamma_weight(v[2], delta2 - 0.25, i, v[2], W) * s)
    end

    delta2 = sum(v2)
    s = line_score(delta2, tl[v[1]], br[v[1]])
    s != 0 && push!(summands, gamma_weight(vsum - 0.5, 2 * delta2 - vsum, v.I..., W) * s)

    for i in v[1]+1:vsum-1
        j = vsum - i
        delta2 = v2[1] + j^2
        s = line_score(delta2, tl[i], br[i])
        s != 0 && push!(summands, gamma_weight(v[1], delta2 - 0.25, v[1], j, W) * s)
    end

    sum_kbn(summands)
end

const Buffer = Array{Int, 3}

mutable struct State
    diags::BitMatrix
    blocked::BitMatrix
    position::Int
    vertices::Vector{Vertex}
    buffer::Buffer
    children::SubArray{Int, 2, Buffer}
    parents::SubArray{Int, 2, Buffer}
    grandparents::SubArray{Int, 2, Buffer}

    function State(n::Int; blocked=nothing)
        isnothing(blocked) && (blocked = falses(n, n))
        checksquare(blocked) == n ||
            throw(DimensionMismatch("size(blocked) $(size(blocked)) != ($n, $n)"))
        buffer = Buffer(undef, 2n+1, n+2, 3)
        children = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
        parents = view(buffer, 1:1, 1:1, 1) .= 0
        new(falses(n, n), blocked, 0, sizehint!(Vertex[], n), buffer, children, parents)
    end
end

Base.isempty(s::State)::Bool = s.position == size(s.buffer, 1) - 2

function Base.popfirst!(s::State)
    @assert !isempty(s)

    n = size(s.buffer, 2) - 2
    i = s.position += 1

    s.grandparents = i > n ? view(s.parents, :, 2:2n+1-i) : s.parents
    s.parents = s.children
    s.children = children = view(s.buffer, 1:i+2, 1:(i < n ? i+2 : 2n-i), mod1(i+2, 3))

    width = size(s.grandparents, 2)
    if i < n
        children[:, 1] .= 0:size(children, 1)-1
        children[:, width+2] .= size(children, 1)-1:-1:0
        children = view(children, :, 2:width+1)
        v = Vertex(i, 1)
    else
        v = Vertex(n, i-n+1)
    end
    resize!(s.vertices, width) .= (v + Vertex(-1, 1) * j for j in 0:width-1)

    (s.grandparents, s.parents, children)
end

function propagate_distances!(i::Int, parents, children)
    depth = size(parents, 1) - 1
    l = view(parents, :, i)
    t = view(parents, :, i+1)
    br = view(children, :, i)
    br[1] = l[1] + 1
    br[depth+2] = t[depth+1] + 1
    view(br, 2:depth+1) .= min.(view(l, 2:depth+1), view(t, 1:depth)) .+ 1
end

function score!(s::State, grandparents, parents, children; W=nothing)::Vector{Float64}
    scores = Vector{Float64}(undef, length(s.vertices))
    Threads.@threads for i in eachindex(s.vertices)
        tl = view(grandparents, :, i)
        br = propagate_distances!(i, parents, children)
        scores[i] = s.blocked[s.vertices[i]] ? -1 : gamma_score(s.vertices[i], tl, br, W)
    end
    scores
end

function setdiags!(s::State, diag_indices, grandparents, children)
    s.diags[s.vertices[diag_indices]] .= true
    Threads.@threads for i in diag_indices
        children[2:s.position+1, i] .= view(grandparents, 1:s.position, i) .+ 1
    end
end

function step!(s::State, (grandparents, parents, children); W=nothing)
    scores = score!(s, grandparents, parents, children, W=W)
    setdiags!(s, findall(>=(0.0), scores), grandparents, children)
end

function grow_gamma_diags(n::Int; W=nothing, margin=0, blocked=nothing)::BitMatrix
    m = n + margin
    s = State(m, blocked=blocked)
    for _ in 1:2m-1
        step!(s, popfirst!(s), W=W)
    end
    s.diags[(margin+1)*onexy:m*onexy]
end

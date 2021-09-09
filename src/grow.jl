const Atom = CartesianIndex{2}
const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

score(de, r) = abs(sqrt(de) - r)

function cost(c::AbstractVector{Int}, a::Atom)
    @assert length(c) == a[1] + a[2] + 1 "$a $(length(c))"
    s = 0.0
    #de = a[2]^2
    for i in 1:a[1]
        #de += 2*i - 1
        de = i^2 + a[2]^2
        s += score(de, c[i+1])
    end
    for i in a[1]+2:a[1]+a[2]
        de = a[1]^2 + (a[1]+a[2]-i+1)^2
        s += score(de, c[i])
    end
    s
end

function blank(l::AbstractVector{Int}, t::AbstractVector{Int})::Vector{Int}
    n = length(l)
    @assert length(t) == n
    regions = Vector{Int}(undef, n + 1)
    regions[1] = l[1] + 1
    regions[2:n] .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1
    regions[end] = t[n] + 1
    regions
end

function diag(tl::AbstractVector{Int}, l::AbstractVector{Int}, t::AbstractVector{Int})
    n = length(tl) + 1
    @assert length(l) == length(t) == n
    regions = Vector{Int}(undef, n + 1)
    regions[1] = l[1] + 1
    regions[2:n] .= view(tl, 1:n-1) .+ 1
    regions[end] = t[n] + 1
    regions
end

function grow(tl, l, t, a)
    bcell = blank(l, t)
    dcell = diag(tl, l, t)
    if cost(dcell, a) < cost(bcell, a)
        (true, dcell)
    else
        (false, bcell)
    end
end

function grow(n)
    diags = Matrix{Union{Missing, Bool}}(missing, n, n)
    grandparents = Dict(Atom(1, 1)=>[0])
    parents = Dict{Atom, Vector{Int}}()
    for (i, r) in enumerate(eg.ribbons(eg.Atom.([(1, 1)]), eg.Atom.([(0, 1), (1, 0)]), n))
        if i <= n
            parents[Atom(i+1, 1)] = collect(0:i)
            parents[Atom(1, i+1)] = collect(i:-1:0)
        end
        children = Dict{Atom, Vector{Int}}()
        for a in r
            diags[a], children[a + onexy] = grow(
                grandparents[a], parents[a + onex], parents[a + oney], a)
        end
        grandparents = parents
        parents = children
    end
    diags
end

function grow2(tl, l, t, a, diags, children, j)
    n = length(l)
    bcell = view(children, :, j)
    bcell[1] = l[1] + 1
    bcell[2:n] .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1
    bcell[end] = t[n] + 1
    dcell = diag(tl, l, t)
    if cost(dcell, a) < cost(bcell, a)
        diags[a] = true
        bcell .= dcell
    else
        diags[a] = false
    end
end


function grow2(n)
    diags = BitMatrix(undef, n, n)
    grandparents = zeros(Int, 1, 1)
    parents = [0 1; 1 0]
    for i in 1:n
        children = Matrix{Int}(undef, i+2, i+2)
        children[:, 1] .= 0:i+1
        a = Atom(i, 1)
        for j in 2:i+1
            tl = view(grandparents, :, j-1)
            l = view(parents, :, j-1)
            t = view(parents, :, j)
            grow2(tl, l, t, a, diags, children, j)
            a += CartesianIndex(-1, 1)
        end
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)
    for i in n-1:-1:1
        children = Matrix{Int}(undef, 2*n-i+2, i)
        a = Atom(n, n-i+1)
        for j in 1:i
            tl = view(grandparents, :, j)
            l = view(parents, :, j)
            t = view(parents, :, j+1)
            grow2(tl, l, t, a, diags, children, j)
            a += CartesianIndex(-1, 1)
        end
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end
    diags
end

function score(g)
    g = BitMatrix(g)
    n = size(g)[1]
    step = Int(n / 8)
    sum(eg.arc(g, step) .!= eg.earc(n, step)) / n^2
end

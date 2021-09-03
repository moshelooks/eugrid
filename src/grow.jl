const Atom = CartesianIndex{2}
const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

score(de, r) = abs(sqrt(de) - r)

function cost(c::Vector{Int}, a::Atom)
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

function blank(l::Vector{Int}, t::Vector{Int})::Vector{Int}
    n = length(l)
    @assert length(t) == n
    regions = Vector{Int}(undef, n + 1)
    regions[1] = l[1] + 1
    regions[2:n] .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1
    regions[end] = t[n] + 1
    regions
end

function diag(tl::Vector{Int}, l::Vector{Int}, t::Vector{Int})
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
    for i in 1:n*2-1
        if i <= n
            parents[Atom(i+1, 1)] = collect(0:i)
            parents[Atom(1, i+1)] = collect(i:-1:0)
        end
        children = Dict{Atom, Vector{Int}}()
        for (a, tl) in grandparents
            maximum(a.I) > n && continue
            diags[a], children[a + onexy] = grow(
                tl, parents[a + onex], parents[a + oney], a)
        end
        grandparents = parents
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

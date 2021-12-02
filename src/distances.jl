const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

const Atom = CartesianIndex{2}

struct DistanceMatrix
    data::Matrix{Int}

    function DistanceMatrix(a::Atom)
        n, m = a.I
        data = Matrix{Int}(undef, n, m)
        data[1:n, m] = n-1:-1:0
        data[n, 1:m-1] = m-1:-1:1
        new(data)
    end
end

atom(d::DistanceMatrix)::Atom = Atom(size(d.data))

function diag(tl::DistanceMatrix)::DistanceMatrix
    br = DistanceMatrix(atom(tl) + onexy)
    br.data[CartesianIndices(tl.data)] .= tl .+ 1
    br
end

function nodiag(t::DistanceMatrix, l::DistanceMatrix)::DistanceMatrix
    assert atom(t) + onex == atom(l) + oney "top $size(t) incompatible with left $size(l)"
    br = DistanceMatrix(atom(t) + onex)
    ix = CartesianIndices(atom(t) - oney)
    br.data[ix] .= min.(view(t.data, ix), view(l.data, ix)) .+ 1
    br
end

isdiag(d::DistanceMatrix)::Bool = 2 - d.data[atom(d) - onexy]

const GraphDistances = Dict{Atom, DistanceMatrix}

struct Extendable
    frozen::GraphDistances
    dcache::IdDict{DistanceMatrix, DistanceMatrix}
    ndcache::Dict{NTuple{2, UInt64}, DistanceMatrix}
    alist::LinkedList{Tuple{Atom, DistanceMatrix}}
end

Extendable(ds::GraphDistances) =
    Extendable(
        ds,
        IdDict{DistanceMatrix, DistanceMatrix}(),
        Dict{NTuple{2, UInt64}, DistanceMatrix}(),
        nil(Tuple{Atom, DistanceMatrix}))

Extendable(n::Int) =
    Extendable(Dict(i=>DistanceMatrix(i)
                    for i in Iterators.flatten(CartesianIndices.(((n, 1), (1, n))))))

graph_distances(e::Extendable) = merge!(Dict(e.alist), e.frozen)

lookup(e::Extendable, a::Atom)::DistanceMatrix =
    get(e.frozen, a) do
        for (i, d) in e.alist
            i == a && return d
        end
    end

function extend(e::Extendable, a::Atom, diag::Bool)::Tuple{Extendable, DistanceMatrix}
    if diag
        tl = e.frozen[a - onexy]
        br = get!(e.dcache, tl) do; diag(tl) end
    else
        t = lookup(e, a - oney)
        l = lookup(e, a - onex)
        br = get!(e.ndcache, (objectid(du), objectid(dl))) do; nodiag(t, l) end
    end
    Extendable(e.frozen, e.dtcache, e.dfcache, cons((a, da), e.alist))
end

function eugrid(e::Extendable)
    diags = Matrix{Union{Missing, Bool}}(missing, maximum(keys(e.frozen)).I)
    for (i, br) in graph_distances(e)
        minimum(i.I) > 1 && (diags[i - onexy] = isdiag(br))
    end
    diags
end

Solver

struct Onion
    layers::Vector{Extendable}

    Onion(n::Int) = new([Extendable(n)])
end

wrap!(o::Onion) = push!(o.layers, Extendable(freeze(o.layers[end])))

peel!(o::Onion) = pop!(o.layers)

eugrid(o::Onion)::BitMatrix = BitMatrix(eugrid(o.layers[end]))

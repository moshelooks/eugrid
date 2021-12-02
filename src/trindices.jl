struct TrilIndices
    n::Int
end

Base.iterate(t::TrilIndices, (i, j)=(1, 1)) =
    j > t.n ? nothing : (CartesianIndex(i, j), i+j>t.n ? (1, j+1) : (i+1, j))

Base.eltype(::Type{TrilIndices}) = CartesianIndex{2}

Base.length(t::TrilIndices) = div(t.n * (t.n + 1), 2)

Base.iterate(rt::Iterators.Reverse{TrilIndices}, (i, j)=(1, rt.itr.n)) =
    j < 1 ? nothing : (CartesianIndex(i, j), i == 1 ? (rt.itr.n + 2 - j, j -1) : (i - 1, j))

struct TrilIndices
    n::Int
end

Base.iterate(t::TriuIndices, (i, j)=(t.n, 1)) =
    j > t.n ? nothing : (CartesianIndex(i, j), i >= t.n ? (1, j+1) : (i+1, j))

Base.eltype(::Type{TriuIndices}) = CartesianIndex{2}

Base.length(t::TriuIndices) = div(t.n * (t.n + 1), 2)

Base.iterate(rt::Iterators.Reverse{TriuIndices}, (i, j)=(1, rt.itr.n)) =
    j < 1 ? nothing : (CartesianIndex(i, j), i == 1 ? (rt.itr.n + 2 - j, j -1) : (i - 1, j))

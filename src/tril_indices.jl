struct TrilIndices2
    n::Int
end

Base.iterate(t::TrilIndices2, (i, j)=(1, 1)) =
    j > t.n ? nothing : (CartesianIndex(i, j), i+j>t.n ? (1, j+1) : (i+1, j))

Base.eltype(::Type{TrilIndices2}) = CartesianIndex{2}

Base.length(t::TrilIndices2) = div(t.n * (t.n + 1), 2)

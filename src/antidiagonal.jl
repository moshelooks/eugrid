antidiagind(m::Integer, n::Integer, k::Integer=0) =
    k <= 0 ?
    range(m+k, step=max(m-1, 1),  length=min(m+k, n)) :
    range(m*(k+1), step=max(m-1, 1), length=min(m, n-k))

function antidiagind(A::AbstractMatrix, k::Integer=0)
    Base.require_one_based_indexing(A)
    antidiagind(size(A)..., k)
end

antidiag(A::AbstractMatrix, k::Integer=0) = A[antidiagind(A, k)]

struct UpperAntiTriangularIndices
    n::Int
end

UpperAntiTriangularIndices(A::AbstractMatrix) =
    UpperAntiTriangularIndices(LinearAlgebra.checksquare(A))

Base.iterate(I::UpperAntiTriangularIndices, (i, j)=(1, 1)) =
    j > I.n ? nothing : (CartesianIndex(i, j), i+j>I.n ? (1, j+1) : (i+1, j))

Base.eltype(::Type{UpperAntiTriangularIndices}) = CartesianIndex{2}

Base.length(I::UpperAntiTriangularIndices) = div(I.n * (I.n + 1), 2)

Base.iterate(r::Iterators.Reverse{UpperAntiTriangularIndices}, (i, j)=(1, r.itr.n)) =
    j < 1 ? nothing : (CartesianIndex(i, j), i == 1 ? (r.itr.n + 2 - j, j -1) : (i - 1, j))

struct LowerAntiTriangularIndices
    n::Int
end

Base.iterate(I::LowerAntiTriangularIndices, (i, j)=(I.n, 1)) =
    j > I.n ? nothing : (CartesianIndex(i, j), i>=I.n ? (I.n-j, j+1) : (i+1, j))

Base.eltype(::Type{LowerAntiTriangularIndices}) = CartesianIndex{2}

Base.length(I::LowerAntiTriangularIndices) = div(I.n * (I.n + 1), 2)

Base.iterate(r::Iterators.Reverse{LowerAntiTriangularIndices}, (i, j)=(r.itr.n, r.itr.n)) =
    j < 1 ? nothing : (CartesianIndex(i, j),
                       i == r.itr.n - j + 1 ? (r.itr.n, j -1) : (i - 1, j))

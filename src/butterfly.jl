struct Butterfly{T, A <: AbstractMatrix{T}}
    data::A

    function Butterfly(data::AbstractMatrix{T}) where {T}
        LinearAlgebra.checksquare(data)
        new{T, typeof(data)}(data)
    end

end

Butterfly{T}(initializer, n::Integer) where {T} = Butterfly(Matrix{T}(initializer, n, n))



struct Spiral
    n::Int
    o::CartesianIndex
    clockwise::Bool
end

const directions = CartesianIndex.([(1, 0), (0, 1), (-1, 0), (0, 1)])

function Base.iterate(s::Spiral, state=())
    at, u, n, m = state
    if n % m == 0
        dix += s.clockwise ? -1 : 1
    end
    if n == 2 * m







    u::CartesianIndex
    n::Int
    m::Int
end

Spiral() = Spiral(CartesianIndex(0, 0), CartesianIndex(1, 0),

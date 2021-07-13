module Eugrid

export Pythagorean, Torus

using CircularArrays: CircularArray

const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

include("Pythagorean.jl")
include("torus.jl")
include("exnil.jl")

end # module

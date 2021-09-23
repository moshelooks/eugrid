module Eugrid

#export Onion, Triple, all_triples, check_diags

import Statistics

using CircularArrays: CircularArray
using LinearAlgebra: checksquare
using StableRNGs: StableRNG

include("geometry.jl")
include("search.jl")
include("validation.jl")
include("grow.jl")
include("analysis.jl")

end # module

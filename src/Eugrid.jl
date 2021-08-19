module Eugrid

using Base.Iterators: flatten

#export Onion, Triple, all_triples, check_diags

#const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

#const Atom = CartesianIndex{2}
#Atom(n::Int) = Atom(n, n)

include("geometry.jl")
#include("constraints.jl")
#include("onion.jl")

end # module

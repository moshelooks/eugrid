module Eugrid

export Grid, Vertex, sps, chessboard, manhattan, vertices, isplanar, eccentricity, geodesics,
    circle_points, midpoints, two_circle_points
export grow_corner_diags, gamma_score, sparsity_cutoff, grow_grid


import DataFrames, GLM, Statistics

using LinearAlgebra: checksquare
using StableRNGs: StableRNG

include("geometry.jl")
include("growth.jl")
#=include("search.jl")
include("validation.jl")

include("analysis.jl")=#

end # module

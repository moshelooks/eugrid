module Eugrid

export Grid, chessboard, manhattan, Vertex, vertices, isplanar, sps, distance, eccentricity,
    euclidean_eccentricity, geodesics, circle_points, midpoints, two_circle_points
export grow_corner_diags, gamma_score, sparsity_cutoff, grow_gamma_diags, grow_grid
export score, sample

import DataFrames, GLM, Statistics

using KahanSummation: sum_kbn
using LinearAlgebra: checksquare
using StableRNGs: StableRNG

include("geometry.jl")
include("growth.jl")
include("sampling.jl")

end # module

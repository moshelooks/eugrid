module Eugrid

export sps, disorder
export b2d, d2b, i2d, d2i
export Grid, chessboard, manhattan, isplanar
export Vertex, vertices
export distance, eccentricity, euclidean_eccentricity, geodesics
export circle_points, midpoints, two_circle_points
export crisscross

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

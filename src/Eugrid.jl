module Eugrid

export sps, disorder
export b2d, d2b, i2d, d2i
export Grid, chessboard, manhattan, randgrid, isplanar
export Vertex, vertices
export distance, eccentricity, euclidean_eccentricity, expected_euclidean_eccentricity
export geodesics, circle_points, midpoints, two_circle_points
export crisscross

export grow_corner_diags, gamma_score, grow_gamma_diags, checkerboard, grow_grid

export score, sample

import DataFrames, GLM, OffsetArrays, Optim, Plots, Statistics

using KahanSummation: sum_kbn
using LinearAlgebra: checksquare
using StableRNGs: StableRNG

include("geometry.jl")
include("growth.jl")
include("sampling.jl")
include("figures.jl")

end # module

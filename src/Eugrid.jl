module Eugrid

# geometry.jl
export sps, disorder
export b2d, d2b, i2d, d2i
export Grid, chessboard, manhattan, isplanar
export Vertex, vertices
export distance, eccentricity, euclidean_eccentricity
export geodesics, circle_points, midpoints, two_circle_points

# growth.jl
export grow_corner_diags, gamma_score, grow_gamma_diags

# sampling.jl
export geodesic_model, diags_sampling, randgrid, rand_sampling, disordered_NG_sampling

# figures.jl
export draw_hg, draw_gg, draw_tg, draw_ag, draw_rand, draw_gamma, draw_gamma_disordered

import DataFrames, GLM, OffsetArrays, Optim, Plots, Statistics

using KahanSummation: sum_kbn
using LaTeXStrings
using LinearAlgebra: checksquare
using StableRNGs: StableRNG

include("geometry.jl")
include("growth.jl")
include("sampling.jl")
include("figures.jl")

end # module

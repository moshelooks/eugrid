using eugrid
using Test

import LinearAlgebra as la

@testset "shortest_paths_symmetric" begin
    diags = Bool.([1 0 0 0
                   0 0 1 0
                   0 1 1 0
                   0 0 0 0])
    d, b = shortest_paths(diags)

    @test d == [0 1 2 3 4
                1 1 2 3 4
                2 2 3 3 4
                3 3 3 4 5
                4 4 4 5 6]

    @test sort.(broadcast.(x->x.I, collect.(b))) == hcat(
        [[], [], [], [], []],
        [[], [(1, 1)], [(1, 1)], [(1, 1)], [(1, 1)]],
        [[], [(1, 1)], [(1, 1)], [(1, 1), (3, 2)], [(1, 1), (3, 2)]],
        [[], [(1, 1)], [(1, 1), (2, 3)], [(1, 1)], [(1, 1)]],
        [[], [(1, 1)], [(1, 1), (2, 3)], [(1, 1)], [(1, 1)]])

end

@testset "shortest_paths_off_diag" begin
    diags = Bool.([0 0 0 0
                   1 0 0 0
                   0 1 0 0
                   0 0 1 0])
    d, b = shortest_paths(diags)

    @test d == [0 1 2 3 4
                1 2 3 4 5
                2 2 3 4 5
                3 3 3 4 5
                4 4 4 4 5]

    @test sort.(broadcast.(x->x.I, collect.(b))) == hcat(
        [[], [], [], [], []],
        [[], [], [(2, 1)], [(2, 1)], [(2, 1)]],
        [[], [], [(2, 1)], [(2, 1), (3, 2)], [(2, 1), (3, 2)]],
        [[], [], [(2, 1)], [(2, 1), (3, 2)], [(2, 1), (3, 2), (4, 3)]],
        [[], [], [(2, 1)], [(2, 1), (3, 2)], [(2, 1), (3, 2), (4, 3)]])

end

@testset "grid_add" begin
    g = Grid(5)
    for ij in CartesianIndices(g.diags)
        eugrid.add(g, ij.I...)
    end

    for (i, l) in enumerate(eugrid.LIndices(g.n_inner))
        for (j, t) in enumerate(eugrid.TIndices(g.n_inner))
            @test g.dg[i, j] == la.norm(l .- t, l[1] > t[1] && l[2] < t[2] ? Inf : 1)
        end
    end
end

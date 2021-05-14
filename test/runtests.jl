using eugrid
using Test

using LinearAlgebra

@testset "shortest_paths_symmetric" begin
    diags = Bool.([
        1 0 0 0
        0 0 1 0
        0 1 1 0
        0 0 0 0
    ])
    d, b = shortest_paths(diags)

    @test d == [
        0 1 2 3 4
        1 1 2 3 4
        2 2 3 3 4
        3 3 3 4 5
        4 4 4 5 6
    ]

    @test sort.(broadcast.(x -> x.I, collect.(b))) == hcat(
        [[], [], [], [], []],
        [[], [(1, 1)], [(1, 1)], [(1, 1)], [(1, 1)]],
        [[], [(1, 1)], [(1, 1)], [(1, 1), (3, 2)], [(1, 1), (3, 2)]],
        [[], [(1, 1)], [(1, 1), (2, 3)], [(1, 1)], [(1, 1)]],
        [[], [(1, 1)], [(1, 1), (2, 3)], [(1, 1)], [(1, 1)]],
    )

end

@testset "shortest_paths_off_diag" begin
    diags = Bool.([
        0 0 0 0
        1 0 0 0
        0 1 0 0
        0 0 1 0
    ])
    d, b = shortest_paths(diags)

    @test d == [
        0 1 2 3 4
        1 2 3 4 5
        2 2 3 4 5
        3 3 3 4 5
        4 4 4 4 5
    ]

    @test sort.(broadcast.(x -> x.I, collect.(b))) == hcat(
        [[], [], [], [], []],
        [[], [], [(2, 1)], [(2, 1)], [(2, 1)]],
        [[], [], [(2, 1)], [(2, 1), (3, 2)], [(2, 1), (3, 2)]],
        [[], [], [(2, 1)], [(2, 1), (3, 2)], [(2, 1), (3, 2), (4, 3)]],
        [[], [], [(2, 1)], [(2, 1), (3, 2)], [(2, 1), (3, 2), (4, 3)]],
    )

end

@testset "helpers" begin
    @test collect(eugrid.lindices(2)) == [(1, 0), (2, 0), (2, 1)]
    @test collect(eugrid.tindices(2)) == [(0, 1), (0, 2), (1, 2)]

    @test eugrid.taxicab([0, 1], [3, 5]) == eugrid.taxicab([3, 5], [0, 1]) == 7
    @test eugrid.euclid([0, 1], [1, 2]) == eugrid.euclid([1, 2], [0, 1]) == sqrt(2)
end

@testset "grid_attrs" begin
    g = Grid(4)
    @test g.n == 4
    @test g.diags == falses(4, 4)
    @test size(g.dil) == size(g.dit) == (4, 4, 7)

    for x = 1:4
        for y = 1:4
            for (i, l) in enumerate(eugrid.lindices(4))
                @test g.dil[x, y, i] == eugrid.taxicab((x, y - 1), l)
            end
            for (i, t) in enumerate(eugrid.tindices(4))
                @test g.dit[x, y, i] == eugrid.taxicab((x - 1, y), t)
            end
        end
    end

    @test g.dg == [
        2 3 4 5 4 5 6
        3 4 5 6 5 4 5
        4 5 6 7 6 5 4
        5 6 7 8 7 6 5
        4 5 6 7 6 5 4
        5 4 5 6 5 4 3
        6 5 4 5 4 3 2
    ]

    @test isapprox(
        g.de .^ 2,
        [
            2 5 10 17 16 17 20
            5 8 13 20 17 16 17
            10 13 18 25 20 17 16
            17 20 25 32 25 20 17
            16 17 20 25 18 13 10
            17 16 17 20 13 8 5
            20 17 16 17 10 5 2
        ],
    )

    @test size(g.dd) == (7, 7)
    @test issymmetric(@. Int(floor(g.dd * 1e9)))
    @test isapprox(g.dd[1, 1], 2 .* (g.dg[1, 1] .- g.de[1, 1]) .- 1)
end

@testset "grid_add_diag" begin
    g = Grid(4)
    for ij in CartesianIndices(g.diags)
        eugrid.add_diag(g, ij.I...)
    end

    expected_d(v, w) = norm(v .- w, (v[1] > w[1]) == (v[2] < w[2]) ? Inf : 1)

    for x = 1:4
        for y = 1:4
            for (i, l) in enumerate(eugrid.lindices(4))
                @test g.dil[x, y, i] == expected_d(l, (x, y - 1))
            end
            for (i, t) in enumerate(eugrid.tindices(4))
                @test g.dit[x, y, i] == expected_d(t, (x - 1, y))
            end
        end
    end

    for (i, l) in enumerate(eugrid.lindices(g.n))
        for (j, t) in enumerate(eugrid.tindices(g.n))
            @test g.dg[i, j] == expected_d(l, t)
            @test issymmetric(@. Int(floor(g.dd * 1e9)))
            @test isapprox(g.dd[1, 1], 2 .* (g.dg[1, 1] .- g.de[1, 1]) .- 1)
        end
    end
end

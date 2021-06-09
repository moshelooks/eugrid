using eugrid
using Test

using LinearAlgebra

@testset "torus_ctor" begin
    t = eugrid.Torus(5, 1)
    @test t.n == 5
    @test t.k == 1
    @test t.m == 3
    @test t.w == 5
    @test t.ls == CartesianIndex.([0, 1, 2], [2, 1, 0])
    @test t.sed == [18 20 26; 20 18 20; 26 20 18]
    @test t.dk == fill(2, (3, 5, 5))
    @test t.dm == fill(6, (5, 5, 5))
    @test t.scores == fill(eugrid.score(t, CartesianIndex(1, 1)), (5, 5))
    @test t.diags == falses(5, 5)
end

@testset "wrap" begin
    t = eugrid.Torus(5, 1)
    unwrapped = [4, 5, 1, 2, 3, 4, 5, 1]
    @test eugrid.wrap.([t], CartesianIndices((-1:6, -1:6))) ==
        CartesianIndex.(unwrapped, reshape(unwrapped, 1, :))
end

raw_score(x, y, d) = d * (3 * (x^2 + y^2) - d^2)

@testset "score" begin
    t = eugrid.Torus(5, 1)
    xs = [3, 2, 1, 4, 3, 2, 5, 4, 3]
    ys = [3, 4, 5, 2, 3, 4, 1, 2, 3]
    d0s = xs .+ ys
    d1s = d0s .- 1
    @test t.scores == fill(sum(raw_score.(xs, ys, d1s) .- raw_score.(xs, ys, d0s)), (5, 5))
end

@testset "staircase" begin
    t = eugrid.Torus(5, 1)
    @test collect(eugrid.Staircase(t, CartesianIndex(1, 1))) ==
        CartesianIndex.([1, 2, 3, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 4, 5, 1],
                        [4, 4, 4, 5, 5, 5, 5, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3])
end

@testset "add_diag" begin
    t = eugrid.Torus(5, 1)
    u = CartesianIndex(1, 1)
    eugrid.add_diag(t, u)
    @test t.diags[u] == true
    @test sum(t.diags) == 1

    @test t.scores[u] == 0
end

#=
    xs = [3, 2, 1, 4, 3, 2, 5, 4, 3]
    ys = [3, 4, 5, 2, 3, 4, 1, 2, 3]
    d0s = xs .+ ys
    d1s = d0s .- 1
    @test t.scores == fill(sum(raw_score.(xs, ys, d1s) .- raw_score.(xs, ys, d0s)), (5, 5))


end



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

    for x = 1:4, y = 1:4
        for (i, l) in enumerate(eugrid.lindices(4))
            @test g.dil[x, y, i] == eugrid.taxicab((x, y - 1), l)
        end
        for (i, t) in enumerate(eugrid.tindices(4))
            @test g.dit[x, y, i] == eugrid.taxicab((x - 1, y), t)
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

    for x = 1:4, y = 1:4
        for (i, l) in enumerate(eugrid.lindices(4))
            @test g.dil[x, y, i] == expected_d(l, (x, y - 1))
        end
        for (i, t) in enumerate(eugrid.tindices(4))
            @test g.dit[x, y, i] == expected_d(t, (x - 1, y))
        end
    end

    for (i, l) in enumerate(eugrid.lindices(g.n)), (j, t) in enumerate(eugrid.tindices(g.n))
        @test g.dg[i, j] == expected_d(l, t)
        @test issymmetric(@. Int(floor(g.dd * 1e9)))
        @test isapprox(g.dd[1, 1], 2 .* (g.dg[1, 1] .- g.de[1, 1]) .- 1)
    end
end


function validate_flip(g, f)
    @assert !g.diags[f.x, f.y]
    lt = Set((l, t) for l = 1:2*g.n-1, t = 1:2*g.n-1)
    for i in f.impacts
        for j in f.impacts
            @test (i === j) == eugrid.intersects(i, j)
        end
        for l in i.l, t in i.t
            @test eugrid.improves(g, f, l, t)
            pop!(lt, (l, t))
        end
    end
    for (l, t) in lt
        @test !eugrid.improves(g, f, l, t)
    end
end

@testset "flipper_attrs" begin
    flipper = eugrid.Flipper(4)
    for f in flipper.flips
        validate_flip(flipper.g, f)
    end
end
=#

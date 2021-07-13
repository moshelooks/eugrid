using Eugrid
using Test

@testset "pythagorean" begin
    @test Pythagorean.children([3, 4, 5]) == [[5, 12, 13], [21, 20, 29], [15, 8, 17]]

    @test Pythagorean.primitive_triples(3) == []
    @test Pythagorean.primitive_triples(4) == [(3, 4, 5)]
    @test Pythagorean.primitive_triples(21) == [
        (3, 4, 5), (5, 12, 13), (8, 15, 17), (20, 21, 29)]
    @test Pythagorean.primitive_triples(84) == [
        (3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25),
        (20, 21, 29), (12, 35, 37), (9, 40, 41), (28, 45, 53),
        (11, 60, 61), (16, 63, 65), (33, 56, 65), (48, 55, 73),
        (13, 84, 85), (36, 77, 85), (39, 80, 89), (65, 72, 97)]

    @test Pythagorean.scaled_triples(4) == [(3, 4, 5), (4, 3, 5)]
    @test Pythagorean.scaled_triples(8) == [(6, 8, 10), (8, 6, 10)]
    @test Pythagorean.scaled_triples(10) == [(6, 8, 10), (8, 6, 10)]
    @test Pythagorean.scaled_triples(12) == [
        (9, 12, 15), (12, 9, 15), (5, 12, 13), (12, 5, 13)]
end

@testset "torus_empty" begin
    t = Torus(5, 4)
    @test t.n == 5
    @test t.k == 4
    @test size(t.d) == (t.k, t.k, t.n, t.n)
    for i in 1:t.k, j in 1:t.k
        @test t.d[i, j, :, :] == fill(i + j, (t.n, t.n))
    end
    @test t.diags == falses(t.n, t.n)
    @test size(t.triples) == (t.k, t.k)
    for i in CartesianIndices(t.triples)
        expected = []
        i[1] <= 3 && i[2] <= 4 && push!(expected, (CartesianIndex(3, 4), 5))
        i[1] <= 4 && i[2] <= 3 && push!(expected, (CartesianIndex(4, 3), 5))
        @test t.triples[i] == expected
    end
    @test CartesianIndices(t) == CartesianIndices((t.n, t.n))

    unwrapped = [4, 5, 1, 2, 3, 4, 5, 1]
    @test Eugrid.wrap.([t], CartesianIndices((-1:6, -1:6))) ==
        CartesianIndex.(unwrapped, reshape(unwrapped, 1, :))

    @test Eugrid.deltas(t) == fill(2, (5, 5, 2))
    @test Eugrid.bound(t) == 2

    for u in CartesianIndices(t)
        for v in CartesianIndices((t.k, t.k))
            for w in CartesianIndex(0,0):CartesianIndex(t.k, t.k)
                @test Eugrid.dvia(t, u, v, w) == sum(v.I) + sum(w.I) - 1
            end
        end
    end
end

@testset "torus_add_diag_one" begin
    for u in CartesianIndices((5, 5))
        t = Torus(5, 4)
        Eugrid.add_diag(t, u)
        @test t.diags[u] == true
        @test sum(t.diags) == 1

        expected = Torus(5,4).d
        for i in CartesianIndices((t.k, t.k))
            v = Eugrid.wrap(t, u - i + Eugrid.onexy)
            for j in i:(Eugrid.onexy * t.k)
                expected[j, v] -= 1
            end
        end
        @test t.d == expected

        @test Eugrid.bound(t) == 2
    end
end

@testset "torus_add_diag_all" begin
    t = Torus(5, 4)
    for u in CartesianIndices(t)
        Eugrid.add_diag(t, u)
    end
    @test t.diags == trues(t.n, t.n)

    for i in 1:t.k, j in 1:t.k
        @test t.d[i, j, :, :] == fill(max(i, j), (t.n, t.n))
    end

    @test Eugrid.deltas(t) == fill(-1, (5, 5, 2))
    @test Eugrid.bound(t) == 1

    for u in CartesianIndices(t)
        for v in CartesianIndices((t.k, t.k))
            for w in CartesianIndex(0,0):CartesianIndex(t.k, t.k)
                @test Eugrid.dvia(t, u, v, w) == max(v.I...) + max(w.I...)
            end
        end
    end
end
#=


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
=#


#=
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
=#

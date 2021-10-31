using Eugrid
using Test

const eg = Eugrid

@testset "grid" begin
    @test_throws DimensionMismatch Grid(falses(2, 3), falses(2, 2))
    @test_throws DimensionMismatch Grid(falses(2, 2), falses(2, 3))
    @test_throws DimensionMismatch Grid(falses(2, 2), falses(3, 3))
    Grid(falses(3, 3), falses(3, 3)).n == 4
end

@testset "isplanar" begin
    @test isplanar(manhattan(4))
    @test isplanar(chessboard(1))
    @test isplanar(manhattan(1))
    @test !isplanar(chessboard(2))
    @test isplanar(manhattan(2))

    diags = falses(4, 4)
    antidiags = falses(4, 4)

    diags[1, 1] = true
    isplanar(Grid(diags, antidiags))

    antidiags[2, 1] = true
    @test isplanar(Grid(diags, antidiags))

    diags[2, 1] = true
    @test !isplanar(Grid(diags, antidiags))
end

@testset "clamp" begin
    g = manhattan(4)
    @test clamp(Vertex(-1, 5), g) == Vertex(1, 4)
    @test clamp(Vertex(5, -1), g) == Vertex(4, 1)
    @test clamp(Vertex(3, 4), g) == Vertex(3, 4)
    @test clamp(Vertex(-2, -3), g) == Vertex(1, 1)
    @test clamp(Vertex(9, 9), g) == Vertex(4, 4)
end

@testset "vertices" begin
    @test vertices(chessboard(4)) == CartesianIndices((4, 4))
end

@testset "sps_distance_chessboard" begin
    g = chessboard(5)
    for v in vertices(g)
        dm = [sps(g, v, m) for m in 1:5]
        for i in vertices(g)
            di = maximum(abs.(i.I .- v.I))
            @test distance(g, v, i) == distance(g, i, v) == i2d(di)
            for (m, d) in enumerate(dm)
                expected = di > m ? typemax(Int32) : di
                @test d[i] == i2d(expected)
            end
        end
    end
end

@testset "sps_distance_manhattan" begin
    g = manhattan(5)
    for v in vertices(g)
        dm = [sps(g, v, m) for m in 1:5]
        for i in vertices(g)
            delta = abs.(i.I .- v.I)
            di = sum(delta)
            @test distance(g, v, i) == distance(g, i, v) == i2d(di)
            for (m, d) in enumerate(dm)
                expected = maximum(delta) > m ? typemax(Int32) : di
                @test d[i] == i2d(expected)
            end
        end
    end
end

@testset "sps_distance_bespoke" begin
    g = Grid(Bool[0 0 1 0;
                  0 1 0 1;
                  0 0 0 0
                  0 0 0 0],
             Bool[0 0 0 0;
                  0 0 1 0;
                  0 1 0 0;
                  0 0 0 0])

    @test sps(g, Vertex(1, 1)) == [distance(g, Vertex(1, 1), v) for v in vertices(g)] ==
        i2d.([0 1 2 3 4;
              1 2 3 3 4;
              2 3 3 4 4;
              3 4 4 5 5;
              4 5 5 6 6])

    @test sps(g, Vertex(3, 3)) == [distance(g, Vertex(3, 3), v) for v in vertices(g)] ==
        i2d.([3 2 2 2 3;
              2 1 1 1 2;
              2 1 0 1 2;
              2 1 1 2 3;
              3 2 2 3 4])

    @test sps(g, Vertex(5, 5)) == [distance(g, Vertex(5, 5), v) for v in vertices(g)] ==
        i2d.([6 5 4 4 4;
              6 5 4 3 3;
              6 5 4 3 2;
              5 4 3 2 1;
              4 3 2 1 0])

    @test sps(g, Vertex(5, 1)) == [distance(g, Vertex(5, 1), v) for v in vertices(g)] ==
        i2d.([4 5 5 5 6
              3 4 4 4 5;
              2 3 3 4 5;
              1 2 3 4 5;
              0 1 2 3 4])

    @test sps(g, Vertex(1, 5)) == [distance(g, Vertex(1, 5), v) for v in vertices(g)] ==
        i2d.([4 3 2 1 0;
              5 4 3 2 1;
              5 4 3 3 2;
              5 4 4 4 3;
              6 5 5 5 4])
end

@testset "eccentricity" begin
    g = Grid(Bool[0 0 1 0;
                  0 1 0 1;
                  0 0 0 0
                  0 0 0 0],
             Bool[0 0 0 0;
                  0 0 1 0;
                  0 1 0 0;
                  0 0 0 0])


    @test eccentricity(g, Vertex(1, 1)) == 6
    @test eccentricity(g, Vertex(3, 3)) == 4
    @test eccentricity(g, Vertex(5, 1)) == 6
    @test eccentricity(g, Vertex(1, 5)) == 6
end

@testset "euclidean_eccentricity" begin
    for n in 2:5, v in vertices(n)
        @test euclidean_eccentricity(n, v) ==
            maximum([sqrt(sum((u-v).I.^2)) for u in vertices(n)])
    end
end

@testset "geodesics" begin
    g = Grid(Bool[0 0 1 0;
                  0 1 0 1;
                  0 0 0 0
                  0 0 0 0],
             Bool[0 0 0 0;
                  0 0 1 0;
                  0 1 0 0;
                  0 0 0 0])

    @test geodesics(g, Vertex(1, 1), Vertex(3, 3)) ==
        geodesics(g, Vertex(3, 3), Vertex(1, 1)) == [1 1 0 0 0;
                                                     1 1 0 0 0;
                                                     0 0 1 0 0;
                                                     0 0 0 0 0;
                                                     0 0 0 0 0]
    @test geodesics(g, Vertex(1, 1), Vertex(5, 5)) ==
        geodesics(g, Vertex(5, 5), Vertex(1, 1)) == [1 1 1 0 0;
                                                     0 0 0 1 0;
                                                     0 0 0 0 1;
                                                     0 0 0 0 1;
                                                     0 0 0 0 1]

    @test geodesics(g, Vertex(5, 1), Vertex(1, 5)) ==
        geodesics(g, Vertex(1, 5), Vertex(5, 1)) == [0 0 0 1 1;
                                                     0 0 0 1 1;
                                                     0 0 1 0 0;
                                                     1 1 0 0 0;
                                                     1 1 0 0 0]


    @test geodesics(g, Vertex(4, 2), Vertex(1, 3)) ==
        geodesics(g, Vertex(1, 3), Vertex(4, 2)) == [0 0 1 0 0;
                                                     0 0 1 1 0;
                                                     0 0 1 0 0;
                                                     0 1 0 0 0;
                                                     0 0 0 0 0]

    @test geodesics(g, Vertex(2, 2), Vertex(3, 5)) ==
        geodesics(g, Vertex(3, 5), Vertex(2, 2)) == [0 0 0 0 0;
                                                     0 1 1 1 0;
                                                     0 0 1 1 1;
                                                     0 0 0 0 0;
                                                     0 0 0 0 0]
end

@testset "circles" begin
    g = Grid(Bool[0 0 1 0;
                  0 1 0 1;
                  0 0 0 0
                  0 0 0 0],
             Bool[0 0 0 0;
                  0 0 1 0;
                  0 1 0 0;
                  0 0 0 0])

    @test circle_points(g, Vertex(1, 1), 0) == Vertex.([(1, 1)])
    @test circle_points(g, Vertex(1, 1), 1) == Vertex.([(2, 1), (1, 2)])
    @test circle_points(g, Vertex(1, 1), 2) == Vertex.([(3, 1), (2, 2), (1, 3)])
    @test circle_points(g, Vertex(1, 1), 3) == Vertex.([
        (4, 1), (3, 2), (2, 3), (3, 3), (1, 4), (2, 4)])


    @test circle_points(g, Vertex(3, 3), 0) == Vertex.([(3, 3)])
    @test circle_points(g, Vertex(3, 3), 1) == Vertex.([
        (2, 2), (3, 2), (4, 2), (2, 3), (4, 3), (2, 4), (3, 4)])
    @test circle_points(g, Vertex(3, 3), 2) == Vertex.([
        (2, 1), (3, 1), (4, 1), (1, 2), (5, 2), (1, 3), (5, 3), (1, 4), (4, 4), (2, 5),
        (3, 5)])

    @test sort(midpoints(g, Vertex(1, 1), Vertex(3, 3))) ==
        sort(midpoints(g, Vertex(3, 3), Vertex(1, 1))) == Vertex.([
            (2, 1), (1, 2), (2, 2)])

    @test sort(midpoints(g, Vertex(1, 1), Vertex(5, 5))) ==
        sort(midpoints(g, Vertex(5, 5), Vertex(1, 1))) == Vertex.([(2, 4)])

    @test two_circle_points(g, Vertex(1, 1), Vertex(3, 3), 1) ==
        two_circle_points(g, Vertex(3, 3), Vertex(1, 1), 1) == Vertex[]

    @test sort(two_circle_points(g, Vertex(1, 1), Vertex(3, 3), 1, strict=false)) ==
        sort(two_circle_points(g, Vertex(3, 3), Vertex(1, 1), 1, strict=false)) == Vertex.([
            (2, 1), (1, 2), (2, 2)])

    @test sort(two_circle_points(g, Vertex(1, 1), Vertex(3, 3), 2)) ==
        sort(two_circle_points(g, Vertex(3, 3), Vertex(1, 1), 2)) == Vertex.([
            (3, 1), (1, 3)])
end

@testset "grow_corner_diags" begin
    n = 11
    le_diags = grow_corner_diags(n-1)
    leq_diags = grow_corner_diags(n-1, true)
    for tl in vertices(n-1)
        d = isqrt(tl[1]^2 + tl[2]^2)

        le_tl = sps(le_diags[eg.onexy:tl-eg.onexy])[end]
        le_l = sps(le_diags[eg.onexy:tl-eg.oney])[end]
        le_t = sps(le_diags[eg.onexy:tl-eg.onex])[end]
        @test le_diags[tl] == (abs(le_tl + 1 - d) < abs(min(le_l, le_t) + 1 - d))

        leq_tl = sps(leq_diags[eg.onexy:tl-eg.onexy])[end]
        leq_l = sps(leq_diags[eg.onexy:tl-eg.oney])[end]
        leq_t = sps(leq_diags[eg.onexy:tl-eg.onex])[end]
        @test leq_diags[tl] == (abs(leq_tl + 1 - d) <= abs(min(leq_l, leq_t) + 1 - d))
    end
end

@testset "gamma_score" begin
    @test gamma_score(Vertex(1, 1), [0], [1]) == 0.0
    @test gamma_score(Vertex(1, 1), [0], [0]) == atan(0.75)
    @test gamma_score(Vertex(1, 1), [2], [0]) == (isqrt(8) - 3) * atan(0.75)

    @test gamma_score(Vertex(2, 1), [0, 0], [1, 1]) == 0.0
    @test gamma_score(Vertex(2, 1), [0, 0], [0, 1]) == atan(4 / 7)
    @test gamma_score(Vertex(2, 1), [2, 0], [0, 1]) == (isqrt(8) - 3) * atan(4 / 7)

    @test gamma_score(Vertex(1, 2), [0, 0], [1, 1]) == 0.0
    @test gamma_score(Vertex(1, 2), [0, 0], [1, 0]) == atan(4 / 7)
    @test gamma_score(Vertex(1, 2), [0, 2], [1, 0]) == (isqrt(8) - 3) * atan(4 / 7)

    @test gamma_score(Vertex(2, 2), [0, 0, 0], [1, 1, 1]) == 0.0

    @test gamma_score(Vertex(3, 1), [0, 0, 0], [1, 1, 1]) == 0.0
    @test gamma_score(Vertex(3, 1), [0, 0, 0], [1, 0, 1]) == 2 * atan(4 / 19)
    @test gamma_score(Vertex(3, 1), [0, 0, 0], [1, 2, 1]) == -2 * atan(4 / 19)
end

@testset "state" begin
    s = eg.State(4)
    @test s.diags == s.blocked == falses(4, 4)
    @test s.position == 0
    @test isempty(s.vertices)
    @test size(s.buffer) == (9, 6, 3)
    @test !isempty(s)

    grandparents, parents, children = popfirst!(s)
    @test s.position == 1
    @test s.vertices == Vertex.([(1, 1)])
    @test grandparents == zeros(Int, 1, 1)
    @test parents == [0 1; 1 0]
    @test size(children) == (3, 1)
    @test size(s.children) == (3, 3)
    @test s.children[:, 1] == [0, 1, 2]
    @test s.children[:, 3] == [2, 1, 0]
    @test !isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [2]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])

    grandparents, parents, children = popfirst!(s)
    @test s.position == 2
    @test s.vertices == Vertex.([(2, 1), (1, 2)])
    @test grandparents == [0 1; 1 0]
    @test parents == [0 1 2; 1 2 1; 2 1 0]
    @test size(children) == (4, 2)
    @test size(s.children) == (4, 4)
    @test s.children[:, 1] == [0, 1, 2, 3]
    @test s.children[:, 4] == [3, 2, 1, 0]
    @test !isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [2, 3]
    @test eg.propagate_distances!(2, parents, children) == [3, 2]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])

    grandparents, parents, children = popfirst!(s)
    @test s.position == 3
    @test s.vertices == Vertex.([(3, 1), (2, 2), (1, 3)])
    @test grandparents == [0 1 2; 1 2 1; 2 1 0]
    @test parents == [0 1 2 3; 1 2 3 2; 2 3 2 1; 3 2 1 0]
    @test size(children) == (5, 3)
    @test size(s.children) == (5, 5)
    @test s.children[:, 1] == [0, 1, 2, 3, 4]
    @test s.children[:, 5] == [4, 3, 2, 1, 0]
    @test !isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [2, 3, 4]
    @test eg.propagate_distances!(2, parents, children) == [3, 4, 3]
    @test eg.propagate_distances!(3, parents, children) == [4, 3, 2]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])

    grandparents, parents, children = popfirst!(s)
    @test s.position == 4
    @test s.vertices == Vertex.([(4, 1), (3, 2), (2, 3), (1, 4)])
    @test grandparents == [0 1 2 3; 1 2 3 2; 2 3 2 1; 3 2 1 0]
    @test parents == [0 1 2 3 4; 1 2 3 4 3; 2 3 4 3 2; 3 4 3 2 1; 4 3 2 1 0]
    @test size(children) == (6, 4)
    @test s.children === children
    @test !isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [2, 3, 4, 5]
    @test eg.propagate_distances!(2, parents, children) == [3, 4, 5, 4]
    @test eg.propagate_distances!(3, parents, children) == [4, 5, 4, 3]
    @test eg.propagate_distances!(4, parents, children) == [5, 4, 3, 2]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])

    grandparents, parents, children = popfirst!(s)
    @test s.position == 5
    @test s.vertices == Vertex.([(4, 2), (3, 3), (2, 4)])
    @test grandparents == [1 2 3; 2 3 4; 3 4 3; 4 3 2; 3 2 1]
    @test parents == [1 2 3 4; 2 3 4 5; 3 4 5 4; 4 5 4 3; 5 4 3 2; 4 3 2 1]
    @test size(children) == (7, 3)
    @test s.children === children
    @test !isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [3, 4, 5, 6, 5]
    @test eg.propagate_distances!(2, parents, children) == [4, 5, 6, 5, 4]
    @test eg.propagate_distances!(3, parents, children) == [5, 6, 5, 4, 3]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])

    grandparents, parents, children = popfirst!(s)
    @test s.position == 6
    @test s.vertices == Vertex.([(4, 3), (3, 4)])
    @test grandparents == [2 3; 3 4; 4 5; 5 4; 4 3; 3 2]
    @test parents == [2 3 4; 3 4 5; 4 5 6; 5 6 5; 6 5 4; 5 4 3; 4 3 2]
    @test size(children) == (8, 2)
    @test s.children === children
    @test !isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [4, 5, 6, 7, 6, 5]
    @test eg.propagate_distances!(2, parents, children) == [5, 6, 7, 6, 5, 4]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])

    grandparents, parents, children = popfirst!(s)
    @test s.position == 7
    @test s.vertices == Vertex.([(4, 4)])
    @test grandparents == reshape([3, 4, 5, 6, 5, 4, 3], :, 1)
    @test parents == [3 4; 4 5; 5 6; 6 7; 7 6; 6 5; 5 4; 4 3]
    @test size(children) == (9, 1)
    @test s.children === children
    @test isempty(s)
    scores = eg.score!(deepcopy((s, grandparents, parents, children))...)
    @test eg.propagate_distances!(1, parents, children) == [5, 6, 7, 8, 7, 6, 5]
    @test isapprox(
        scores, [gamma_score(s.vertices[i], grandparents[:, i], children[2:end-1, i])
                 for i in eachindex(s.vertices)])
end
#=
@testset "step!" begin
    goldens = Bool[1 1 1 1; 1 0 0 0; 1 0 0 1; 1 0 1 1]
    s = eg.State(4)

    for _ in 1:7
        _, _, children = generations = popfirst!(s)
        eg.step!(s, generations)
        for (i, v) in enumerate(s.vertices)
            @test s.diags[v] == goldens[v]
            dv = sps(s.diags[v:-eg.onexy:eg.onexy])
            @test dv[:, end] == children[1:v[1]+1, i]
            @test dv[end, :] == children[end:-1:v[1]+1, i]
        end
    end

    @test s.diags == goldens
end

@testset "grow_gamma_diags" begin
    goldens = Bool[1 1 1 1; 1 0 0 0; 1 0 0 1; 1 0 1 1]
    @test grow_gamma_diags(4, margin=0) == goldens

    @test grow_gamma_diags(4, margin=4) == grow_gamma_diags(8, margin=0)[5:8, 5:8]
end

@testset "grow_grid" begin
    goldens = Bool[1 1 1 1; 1 0 0 0; 1 0 0 1; 1 0 1 1]
    antigoldens = Bool[0 0 0 0; 0 1 0 1; 0 1 1 0; 0 1 0 0]
    g = grow_grid(4, margin=0)
    @test g.diags == goldens
    @test g.antidiags == antigoldens
    @test isplanar(g)

    @test grow_grid(4, margin=4).diags == grow_gamma_diags(8, margin=0)[5:8, 5:8]
    @test grow_grid(4, margin=4).antidiags == grow_grid(8, margin=0).antidiags[5:8, 5:8]

    for sparsity in 0:0.1:1
        gs = grow_grid(4, margin=0, sparsity=sparsity)
        @test all(gs.diags .<= g.diags)
        if sparsity == 0
            @test gs.diags == g.diags
            @test gs.antidiags == g.antidiags
        end
        g = gs
    end
    @test g.diags == g.antidiags == falses(4, 4)
end
=#
#@testset

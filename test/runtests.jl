using Eugrid
using Test

const eg = Eugrid

@testset "grid" begin
    @test_throws DimensionMismatch Grid(falses(2, 3), falses(2, 2))
    @test_throws DimensionMismatch Grid(falses(2, 2), falses(2, 3))
    @test_throws DimensionMismatch Grid(falses(2, 2), falses(3, 3))
    Grid(falses(3, 3), falses(3, 3)).n == 4
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

@testset "isplanar" begin
    g = manhattan(4)
    @test isplanar(g)
    g.diags[1, 1] = true
    @test isplanar(g)
    g.antidiags[2, 1] = true
    @test isplanar(g)
    g.diags[2, 1] = true
    @test !isplanar(g)
    @test isplanar(chessboard(1))
    @test isplanar(manhattan(1))
    @test !isplanar(chessboard(2))
    @test isplanar(manhattan(2))
end

@testset "sps_chessboard" begin
    g = chessboard(5)
    for v in vertices(g), m in 1:5
        d = sps(g, v, m)
        for i in vertices(g)
            di = maximum(abs.(i.I .- v.I))
            expected = di > m ? typemax(Int) : di
            @test d[i] == expected
        end
    end
end

@testset "sps_manhattan" begin
    g = manhattan(5)
    for v in vertices(g), m in 1:5
        d = sps(g, v, m)
        for i in vertices(g)
            di = abs.(i.I .- v.I)
            expected = maximum(di) > m ? typemax(Int) : sum(di)
            @test d[i] == expected
        end
    end
end

@testset "sps_bespoke" begin
    g = Grid([0 0 1 0;
              0 1 0 1;
              0 0 0 0
              0 0 0 0],
             [0 0 0 0;
              0 0 1 0;
              0 1 0 0;
              0 0 0 0])

    @test sps(g, Vertex(1, 1)) == [0 1 2 3 4;
                                   1 2 3 3 4;
                                   2 3 3 4 4;
                                   3 4 4 5 5;
                                   4 5 5 6 6]

    @test sps(g, Vertex(3, 3)) == [3 2 2 2 3;
                                   2 1 1 1 2;
                                   2 1 0 1 2;
                                   2 1 1 2 3;
                                   3 2 2 3 4]

    @test sps(g, Vertex(5, 5)) == [6 5 4 4 4;
                                   6 5 4 3 3;
                                   6 5 4 3 2;
                                   5 4 3 2 1;
                                   4 3 2 1 0]

    @test sps(g, Vertex(5, 1)) == [4 5 5 5 6
                                   3 4 4 4 5;
                                   2 3 3 4 5;
                                   1 2 3 4 5;
                                   0 1 2 3 4]

    @test sps(g, Vertex(1, 5)) == [4 3 2 1 0;
                                   5 4 3 2 1;
                                   5 4 3 3 2;
                                   5 4 4 4 3;
                                   6 5 5 5 4]
end

@testset "eccentricity" begin
    g = Grid([0 0 1 0;
              0 1 0 1;
              0 0 0 0
              0 0 0 0],
             [0 0 0 0;
              0 0 1 0;
              0 1 0 0;
              0 0 0 0])


    @test eccentricity(g, Vertex(1, 1)) == 6
    @test eccentricity(g, Vertex(3, 3)) == 4
    @test eccentricity(g, Vertex(5, 1)) == 6
    @test eccentricity(g, Vertex(1, 5)) == 6
end


@testset "geodesics" begin
    g = Grid([0 0 1 0;
              0 1 0 1;
              0 0 0 0
              0 0 0 0],
             [0 0 0 0;
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
    g = Grid([0 0 1 0;
              0 1 0 1;
              0 0 0 0
              0 0 0 0],
             [0 0 0 0;
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

    @test sort(two_circle_points(g, Vertex(1, 1), Vertex(3, 3), 2)) ==
        sort(two_circle_points(g, Vertex(3, 3), Vertex(1, 1), 2)) == Vertex.([
            (3, 1), (1, 3)])
end

@testset "grow_corner_diags" begin
    n = 11
    le_diags = grow_corner_diags(n-1)
    leq_diags = grow_corner_diags(n-1, true)
    for tl in vertices(n-1)
        d = sqrt(tl[1]^2 + tl[2]^2)

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
    @test gamma_score(Vertex(1, 1), [2], [0]) == (2 * sqrt(2) - 3) * atan(0.75)

    @test gamma_score(Vertex(2, 1), [0, 0], [1, 1]) == 0.0
    @test gamma_score(Vertex(2, 1), [0, 0], [0, 1]) == atan(4 / 7)
    @test gamma_score(Vertex(2, 1), [2, 0], [0, 1]) == (2 * sqrt(2) - 3) * atan(4 / 7)

    @test gamma_score(Vertex(1, 2), [0, 0], [1, 1]) == 0.0
    @test gamma_score(Vertex(1, 2), [0, 0], [1, 0]) == atan(4 / 7)
    @test gamma_score(Vertex(1, 2), [0, 2], [1, 0]) == (2 * sqrt(2) - 3) * atan(4 / 7)

    @test gamma_score(Vertex(2, 2), [0, 0, 0], [1, 1, 1]) == 0.0

    @test gamma_score(Vertex(3, 1), [0, 0, 0], [1, 1, 1]) == 0.0
    @test gamma_score(Vertex(3, 1), [0, 0, 0], [1, 0, 1]) == 2 * atan(4 / 19)
    @test gamma_score(Vertex(3, 1), [0, 0, 0], [1, 2, 1]) == -2 * atan(4 / 19)
end

@testset "sparsity_cutoff" begin
    @test sparsity_cutoff(1:100, 0.0) == 1
    @test sparsity_cutoff(1:100, 0.01) == 2
    @test sparsity_cutoff(1:100, 0.5) == 51
    @test sparsity_cutoff(1:100, 0.98) == 99
    @test sparsity_cutoff(1:100, 0.99) == 100
    @test sparsity_cutoff(1:100, 1.0) > 100
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

    @test grow_gamma_diags(4) == grow_gamma_diags(8, margin=0)[5:8, 5:8]
end


@testset "grow_grid" begin
    goldens = Bool[1 1 1 1; 1 0 0 0; 1 0 0 1; 1 0 1 1]
    antigoldens = Bool[0 0 0 0; 0 1 0 1; 0 1 1 0; 0 1 0 0]
    g = grow_grid(4, margin=0)
    @test g.diags == goldens
    @test g.antidiags == antigoldens
    @test isplanar(g)

    @test grow_grid(4).diags == grow_gamma_diags(8, margin=0)[5:8, 5:8]
    @test grow_grid(4).antidiags == grow_grid(8, margin=0).antidiags[5:8, 5:8]

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



#=

    @test sps(gm, Vertex(1, 1)) == [0 1 2 3 4;
                                    1 2 3 4 5;
                                    2 3 4 4 5;
                                    3 4 4 5 6;
                                    4 5 5 6 7]

    @test eg.DistanceMatrix(eg.Atom(1, 1)).data == zeros(Int, 1, 1)

    d = eg.DistanceMatrix(eg.Atom(2, 3))
    @test d.data[:, 3] == [1, 0]
    @test d.data[2, :] == [2, 1, 0]

    d = eg.DistanceMatrix(eg.Atom(1, 3))
    @test size(d) == (1, 3)
    @test d.a == eg.Atom(1, 3)
    @test eg.Atoms(d) == eg.Atoms(1, 3) == CartesianIndices((1, 3))

    @test eg.diag(d).data == [3 2 1 1; 3 2 1 0]

    l = eg.DistanceMatrix(eg.Atom(2, 2))
    l.data[1, 1] = 42
    @test eg.nodiag(d, l).data == [3 2 1; 2 1 0]
    @test eg.Atoms(l) == eg.Atoms(l.a) == eg.Atoms(2)

    l.data[1, 1] = -3
    @test eg.nodiag(d, l).data == [-2 2 1; 2 1 0]
end

@testset "triple" begin
    @test_throws InexactError eg.Triple(1, 1)

    t = eg.Triple(3, 4)
    @test t.c == 5
    @test t.v == eg.Atom(3, 4)

    t = eg.Triple(6, 8)
    @test t.c == 10
    @test t.v == eg.Atom(6, 8)

    @test eg.all_triples(24) == eg.Triple.(
        [4, 3, 12, 8, 24, 6, 15, 12, 24, 5, 9, 16, 8, 20, 12, 24, 15, 21, 20, 7, 10, 18],
        [3, 4, 5, 6, 7, 8, 8, 9, 10, 12, 12, 12, 15, 15, 16, 18, 20, 20, 21, 24, 24, 24],
    )
end

@testset "region" begin
    t = eg.Triple(3, 4)

    r = eg.Region(t, eg.Atom(1, 1))
    @test r.t === t
    @test r.w == eg.Atom(4, 5)
    @test eg.delta_min(r, eg.Atom(1, 1)) == 4
    @test eg.delta_max(r, eg.Atom(1, 1)) == 7
    @test eg.delta_min(r, eg.Atom(3, 2)) == 3
    @test eg.delta_min(r, eg.Atom(2, 3)) == 2
    @test eg.delta_max(r, eg.Atom(3, 2)) == eg.delta_max(r, eg.Atom(2, 3)) == 4
    @test eg.delta_min(r, eg.Atom(4, 5)) == 0
    @test eg.delta_max(r, eg.Atom(4, 5)) == 0
    for a in eg.Atoms(3, 4)
        @test a in r
    end
    @test !(eg.Atom(4, 1) in r)
    @test !(eg.Atom(1, 5) in r)

    dt = eg.diag(eg.DistanceMatrix(eg.Atom(1, 1)))
    @test eg.delta_distance(r, dt) == 4
    dl = eg.DistanceMatrix(eg.Atom(3, 1))
    @test eg.delta_distance(r, dl) == 3
    d = eg.nodiag(dt, dl)
    @test eg.delta_distance(r, d) == 3

    r = eg.Region(t, eg.Atom(2, 3))
    @test r.t === t
    @test r.w == eg.Atom(5, 7)
    @test eg.delta_min(r, eg.Atom(2, 3)) == 4
    @test eg.delta_max(r, eg.Atom(2, 3)) == 7
    @test eg.delta_min(r, eg.Atom(4, 4)) == 3
    @test eg.delta_min(r, eg.Atom(3, 5)) == 2
    @test eg.delta_max(r, eg.Atom(4, 4)) == eg.delta_max(r, eg.Atom(3, 5)) == 4
    @test eg.delta_min(r, eg.Atom(5, 7)) == 0
    @test eg.delta_max(r, eg.Atom(5, 7)) == 0
    for a in eg.Atoms(2:4, 3:6)
        @test a in r
    end
    @test !(eg.Atom(1, 3) in r)
    @test !(eg.Atom(2, 2) in r)
    @test !(eg.Atom(4, 7) in r)
    @test !(eg.Atom(5, 6) in r)

    dt = eg.DistanceMatrix(eg.Atom(1, 3))
    dl = eg.diag(eg.DistanceMatrix(eg.Atom(1, 1)))
    d = eg.nodiag(dt, dl)
    @test eg.delta_distance(r, d) == 5
end

@testset "regions_of_interest" begin
    b = eg.Box(3)
    @test isempty(b.triples)

    b = eg.Box(3, 4)
    t = eg.Triple(3, 4)
    @test b.triples == [t]
    for a in eg.Atoms(3, 4)
        @test only(eg.regions_of_interest(b, a)) == eg.Region(t, eg.Atom(1, 1))
    end

    roi(b, x, y) = collect(eg.regions_of_interest(b, eg.Atom(x, y)))

    b = eg.Box(4)
    @test b.triples == eg.all_triples(4)
    t43, t34 = b.triples
    @test roi(b, 1, 1) == eg.Region.(b.triples, [eg.Atom(1, 1)])
    @test roi(b, 2, 1) == eg.Region.([t43, t34, t34], eg.Atom.([(1, 1), (1, 1), (2, 1)]))
    @test roi(b, 1, 2) == eg.Region.([t43, t43, t34], eg.Atom.([(1, 1), (1, 2), (1, 1)]))
    @test roi(b, 4, 4) == eg.Region.([t43, t34], eg.Atom.([(1, 2), (2, 1)]))

    b = eg.Box(6, 8, [t])
    roi(x, y) = [r.u for r in eg.regions_of_interest(b.br, t, eg.Atom(x, y))]

    for i in 1:3, j in 1:4
        @test roi(i, j) == eg.Atoms(i, j)
    end

    for i in 1:3, j in 5:8
        @test roi(i, j) == eg.Atoms(i, j-3:5)
    end

    for i in 4:6, j in 1:4
        @test roi(i, j) == eg.Atoms(i-2:4, 1:j)
    end

    for i in 4:6, j in 5:8
        @test roi(i, j) == eg.Atoms(i-2:4, j-3:5)
    end

    b = eg.Box(6, 8, [t, t])
    @test collect(r.u for r in eg.regions_of_interest(b, eg.Atom(3, 2))) ==
        collect(Iterators.flatten((roi(3, 2), roi(3, 2))))

end

@testset "constraints" begin
    t = eg.Triple(3, 4)
    r1 = eg.Region(t, eg.Atom(1, 1))
    r2 = eg.Region(t, eg.Atom(2, 2))
    r3 = eg.Region(t, eg.Atom(3, 3))

    cs = eg.Constraints()
    @test eg.issatisfiable(cs)
    @test isempty(eg.violations(cs))
    @test isempty(eg.clauses(cs))
    @test !eg.isfree(cs, r1)

    eg.cover!(cs, r1)
    eg.cover!(cs, r2)
    eg.free!(cs, r3)
    @test !eg.issatisfiable(cs)
    @test eg.violations(cs) in ([r1, r2], [r2, r1])
    @test collect(eg.clauses(cs)) == [[], []]

    push!(eg.cover!(cs, r1), eg.Atom(4, 4))
    push!(eg.cover!(cs, r2), eg.Atom(5, 5))
    @test eg.issatisfiable(cs)
    @test isempty(eg.violations(cs))
    @test sort(collect(eg.clauses(cs))) == [[eg.Atom(4, 4)], [eg.Atom(5, 5)]]
end
#=
@testset "constrain!" begin
    t = eg.Triple(3, 4)
    b = eg.Box(6, 8, [t])

    cs = eg.Constraints()
    d = eg.diag(eg.DistanceMatrix(eg.Atom(1, 1)))
    eg.constrain!(cs, b, d)
    @test cs.domain == Set([eg.Atom(2, 2)])
    @test cs.region_clauses == Dict(
        eg.Region(t, eg.Atom(1, 1))=>eg._free,
        eg.Region(t, eg.Atom(2, 1))=>eg._free,
        eg.Region(t, eg.Atom(1, 2))=>[eg.Atom(2, 2)],
        eg.Region(t, eg.Atom(2, 2))=>eg._free)
    @test isempty(eg.violations(cs))
    #@test isempty(eg.f(cs))

    cs = eg.Constraints()
    d = eg.diag(d)
    eg.constrain!(cs, b, d)
    @test isempty(cs.domain)
    @test cs.region_clauses == Dict(
        eg.Region(t, eg.Atom(1, 1))=>eg._free,
        eg.Region(t, eg.Atom(2, 1))=>eg._free,
        eg.Region(t, eg.Atom(3, 1))=>[],
        eg.Region(t, eg.Atom(1, 2))=>[],
        eg.Region(t, eg.Atom(2, 2))=>eg._free,
        eg.Region(t, eg.Atom(3, 2))=>eg._free,
        eg.Region(t, eg.Atom(1, 3))=>[],
        eg.Region(t, eg.Atom(2, 3))=>[],
        eg.Region(t, eg.Atom(3, 3))=>eg._free)

    cs = eg.Constraints()
    d = eg.diag(eg.nodiag(eg.DistanceMatrix.([eg.Atom(1, 2), eg.Atom(2, 1)])...))
    eg.constrain!(cs, b, d)
    @test cs.domain == Set([eg.Atom(3, 3)])
    @test cs.region_clauses == Dict(
        eg.Region(t, eg.Atom(1, 1))=>[eg.Atom(3, 3)],
        eg.Region(t, eg.Atom(2, 1))=>eg._free,
        eg.Region(t, eg.Atom(3, 1))=>[eg.Atom(3, 3)],
        eg.Region(t, eg.Atom(1, 2))=>[eg.Atom(3, 3)],
        eg.Region(t, eg.Atom(2, 2))=>eg._free,
        eg.Region(t, eg.Atom(3, 2))=>eg._free,
        eg.Region(t, eg.Atom(1, 3))=>[],
        eg.Region(t, eg.Atom(2, 3))=>[eg.Atom(3, 3)],
        eg.Region(t, eg.Atom(3, 3))=>eg._free)
end

function square(diags, n)
    ds = eg.GraphDistances()
    for i in eg.Atoms(n)
        if minimum(i.I) == 1
            ds[i] = eg.DistanceMatrix(i)
        elseif i - eg.onexy in diags
            ds[i] = eg.diag(ds[i - eg.onexy])
        else
            ds[i] = eg.nodiag(ds[i - eg.onex], ds[i - eg.oney])
        end
    end
    ds
end

@testset "constraints" begin
    t = eg.Triple(3, 4)
    b = eg.Box(6, 8, [t])
    ds = filter(square(Set([eg.Atom(1, 1), eg.Atom(2, 2)]), 3)) do (a, _)
        maximum(a.I) == 3
    end
    cs = eg.constraints(b, ds)
    @test cs.domain == Set(eg.Atom.([(3, 1), (3, 2), (1, 3), (2, 3)]))
    expected = Dict(eg.Region(t, a)=>eg._free for a in eg.Atoms(3))
    @test cs.region_clauses == expected

    ds = filter(square(Set([eg.Atom(2, 2)]), 3)) do (a, _)
        maximum(a.I) == 3
    end
    cs = eg.constraints(b, ds)
    @test cs.domain == Set(eg.Atom.([(3, 1), (3, 2), (1, 3), (2, 3), (3, 3)]))
    expected[eg.Region(t, eg.Atom(1, 1))] = eg.Atom.([(1, 3), (2, 3), (3, 3)])
    foreach(sort!, values(cs.region_clauses))
    @test cs.region_clauses == expected
end

@testset "ribbons" begin
    kernel = eg.Atom.([(1, 1)])
    basis = eg.Atom.([(0, 1), (1, 1), (1, 0)])
    @test eg.ribbons(kernel, basis, 0) == []
    @test eg.ribbons(kernel, basis, 1) == [eg.Atom.([(1, 1)])]
    @test eg.ribbons(kernel, basis, 2) == [
        eg.Atom.([(1, 1)]),
        eg.Atom.([(2, 1), (1, 2), (2, 2)])]
    @test eg.ribbons(kernel, basis, 3) == [
        eg.Atom.([(1, 1)]),
        eg.Atom.([(2, 1), (1, 2), (2, 2)]),
        eg.Atom.([(3, 1), (3, 2), (1, 3), (2, 3), (3, 3)])]

    kernel = eg.Atom.([(1, 2), (2, 1), (1, 1)])
    basis = eg.Atom.([(0, 2), (0, 1), (1, 1), (1, 0), (2, 0)])
    @test eg.ribbons(kernel, basis, 0) == []
    @test eg.ribbons(kernel, basis, 1) == [eg.Atom.([(1, 1)])]
    @test eg.ribbons(kernel, basis, 2) == [
        eg.Atom.([(1, 1), (2, 1), (1, 2)]),
        eg.Atom.([(2, 2)])]
    @test eg.ribbons(kernel, basis, 3) == [
        eg.Atom.([(1, 1), (2, 1), (1, 2)]),
        eg.Atom.([(3, 1), (2, 2), (3, 2), (1, 3), (2, 3)]),
        eg.Atom.([(3, 3)])]
    @test eg.ribbons(kernel, basis, 4) == [
        eg.Atom.([(1, 1), (2, 1), (1, 2)]),
        eg.Atom.([(3, 1), (4, 1), (2, 2), (3, 2), (1, 3), (2, 3), (1, 4)]),
        eg.Atom.([(4, 2), (3, 3), (4, 3), (2, 4), (3, 4)]),
        eg.Atom.([(4, 4)])]
end

@testset "membrane" begin
    exterior = eg.Atom.([(1, 1), (2, 1), (3, 1), (1, 2), (1, 3)])
    m = eg.Membrane(eg.GraphDistances(), exterior)
    @test sort(collect(keys(m.border))) == exterior
    @test isempty(m.targets)

    ds = eg.exterior_distances(m, Set{eg.Atom}())
    @test ds == m.border
    @test isempty(m.diag_cache)
    @test isempty(m.nodiag_cache)

    exterior = eg.Atom.([(2, 2), (3, 2), (2, 3)])
    m = eg.Membrane(ds, exterior)
    @test isempty(m.border)
    @test m.targets == exterior

    diags = Set{eg.Atom}()
    ds1 = eg.exterior_distances(m, diags)
    @test isempty(m.diag_cache)
    @test sort(collect(keys(ds1))) == exterior
    @test ds1[eg.Atom(2, 2)].data == [2 1; 1 0]
    @test ds1[eg.Atom(3, 2)].data == [3 2; 2 1; 1 0]
    @test ds1[eg.Atom(2, 3)].data == [3 2 1; 2 1 0]
    @test sort(objectid.(values(ds1))) == sort(objectid.(values(m.nodiag_cache)))
    @test ds1 == eg.exterior_distances(m, diags)

    diags = Set([eg.Atom(1, 1)])
    ds2 = eg.exterior_distances(m, diags)
    @test only(m.diag_cache) == (eg.Atom(2, 2) => ds2[eg.Atom(2, 2)])
    @test ds2[eg.Atom(2, 2)].data == [1 1; 1 0]
    @test ds2[eg.Atom(3, 2)].data == [2 2; 2 1; 1 0]
    @test ds2[eg.Atom(2, 3)].data == [2 2 1; 2 1 0]
    @test ds2 == eg.exterior_distances(m, diags)

    diags = Set([eg.Atom(2, 1)])
    ds3 = eg.exterior_distances(m, diags)
    @test m.diag_cache[eg.Atom(3, 2)] === ds3[eg.Atom(3, 2)]
    @test ds3[eg.Atom(2, 2)] === ds1[eg.Atom(2, 2)]
    @test ds3[eg.Atom(3, 2)].data == [2 2; 1 1; 1 0]
    @test ds3[eg.Atom(2, 3)]  === ds1[eg.Atom(2, 3)]
    @test ds3 == eg.exterior_distances(m, diags)

    diags = Set([eg.Atom(1, 2)])
    ds4 = eg.exterior_distances(m, diags)
    @test m.diag_cache[eg.Atom(2, 3)] === ds4[eg.Atom(2, 3)]
    @test ds4[eg.Atom(2, 2)] === ds1[eg.Atom(2, 2)]
    @test ds4[eg.Atom(3, 2)] === ds1[eg.Atom(3, 2)]
    @test ds4[eg.Atom(2, 3)].data == [2 1 1; 2 1 0]
    @test ds4 == eg.exterior_distances(m, diags)
end

function validate(a::eg.Assignment)::eg.Assignment
    for c in a.clauses
        for l in c
            @test l in a.free
        end
    end
    for l in a.affirmed
        @test !(l in a.free)
    end
    a
end

function pose(clauses, n)
    atoms = [eg.Atom(i, 1) for i in 1:n]
    cs = eg.Constraints()
    push!(cs.domain, atoms...)
    for (i, c) in enumerate(clauses)
        ix = eg.Region(eg.Triple(3, 4), eg.Atom(i, i))
        cs.region_clauses[ix] = [atoms[j] for j in c]
    end
    validate(eg.Assignment(cs))
end

grovel(atoms) = sort([l[1] for l in atoms])

grovel(a::eg.Assignment) =
    (sort(map(grovel, validate(a).clauses)), grovel(a.free), grovel(a.affirmed))

@testset "assignment" begin
    function simplify(clauses, n)
        a = pose(clauses, n)

        forks = Dict()
        for l in sort(collect(a.free))
            af = deepcopy(a)
            aff = eg.fork!(af, l)
            forks[l[1]] = (grovel(af), grovel(aff))
        end

        (grovel(a)..., forks)
    end

    clauses, free, affirmed, forks = simplify([], 2)
    @test clauses == []
    @test free == [1, 2]
    @test affirmed == []
    @test forks == Dict(
        1=>(([], [2], []), ([], [2], [1])),
        2=>(([], [1], []), ([], [1], [2])))

    clauses, free, affirmed, forks = simplify([[1], [1, 2], [1, 2]], 2)
    @test clauses == []
    @test free == [2]
    @test affirmed == [1]
    @test forks == Dict(
        2=>(([], [], [1]), ([], [], [1, 2])))

    clauses, free, affirmed, forks = simplify([[1, 2], [3], [3], [1, 2]], 3)
    @test clauses == [[1, 2]]
    @test free == [1, 2]
    @test affirmed == [3]
        @test forks == Dict(
        1=>(([], [], [2, 3]), ([], [2], [1, 3])),
        2=>(([], [], [1, 3]), ([], [1], [2, 3])))

    clauses, free, affirmed, forks = simplify([[1, 2], [1, 2, 3], [1, 4]], 4)
    @test clauses == [[1, 2], [1, 4]]
    @test free == [1, 2, 3, 4]
    @test affirmed == []
    @test forks == Dict(
        1=>(([], [3], [2, 4]), ([], [2, 3, 4], [1])),
        2=>(([], [3, 4], [1]), ([[1, 4]], [1, 3, 4], [2])),
        3=>(([[1, 2], [1, 4]], [1, 2, 4], []), ([[1, 2], [1, 4]], [1, 2, 4], [3])),
        4=>(([], [2, 3], [1]), ([[1, 2]], [1, 2, 3], [4])))

    clauses, free, affirmed, forks = simplify([[1, 2, 3], [1, 2, 4], [3]], 4)
    @test clauses == [[1, 2, 4]]
    @test free == [1, 2, 4]
    @test affirmed == [3]
    @test forks == Dict(
        1=>(([[2, 4]], [2, 4], [3]), ([], [2, 4], [1, 3])),
        2=>(([[1, 4]], [1, 4], [3]), ([], [1, 4], [2, 3])),
        4=>(([[1, 2]], [1, 2], [3]), ([], [1, 2], [3, 4])))
end

function solutions(s::eg.Solver)
    xs = [sort(collect(x)) for x in eg.all_solutions!(s)]
    @test length(xs) == length(unique(xs))
    sort!(xs)
end

function solutions(clauses, n, expected...)
    xs = map(grovel, solutions(eg.Solver([pose(clauses, n)])))
    @test xs == sort([sort(Vector{Int}(c)) for c in expected])
end

@testset "solver" begin
    solutions([], 1, [], [1])

    solutions([], 2, [], [1], [2], [1, 2])
    solutions([[1]], 1, [1])
    solutions([[1]], 2, [1], [1, 2])
    solutions([[1], [2]], 2, [1, 2])
    solutions([[1, 2]], 2, [1], [2], [1, 2])
    solutions([[1, 2]], 4, [1], [2], [1, 2],
              [1, 3], [2, 3], [1, 2, 3],
              [1, 4], [2, 4], [1, 2, 4],
              [1, 3, 4], [2, 3, 4], [1, 2, 3, 4])
end

@testset "onion_ctor" begin
    kernel = eg.Atom.([(1, 2), (2, 1), (1, 1)])
    basis = eg.Atom.([(0, 2), (0, 1), (1, 1), (1, 0), (2, 0)])
    o = eg.Onion(kernel, basis, 4)
    @test o.ribbons == eg.ribbons(kernel, basis, 4)
    @test length(o.membranes) == length(o.solvers) == 1
    @test sort(collect(keys(o.membranes[1].interior))) == eg.Atom.([(1, 1), (2, 1), (1, 2)])
    @test sort(collect(keys(o.membranes[1].border))) ==
        eg.Atom.([(3, 1), (4, 1), (1, 3), (1, 4)])
    @test o.membranes[1].targets == eg.Atom.([(2, 2), (3, 2), (2, 3)])
    @test length(o.solvers[1].stack) == 1
    @test isempty(o.solvers[1].stack[1].clauses)
    @test o.solvers[1].stack[1].free == Set(eg.Atom.([(1, 1), (2, 1), (1, 2)]))
    @test isempty(o.solvers[1].stack[1].affirmed)
    @test isempty(o.diags)
    @test o.box.triples == eg.Triple.([4, 3], [3, 4])
end

@testset "onion_step!" begin
    kernel = eg.Atom.([(1, 1)])
    basis = eg.Atom.([(0, 1), (1, 1), (1, 0)])
    o = eg.Onion(kernel, basis, 2)
    @test isequal(eg.eugrid(o), [missing missing; missing missing])

    eg.step!(o)
    @test isequal(eg.eugrid(o), [false missing; missing missing])

    eugrids = map(1:8) do _
        eg.step!(o)
        eg.eugrid(o)
    end
    @test Set(eugrids) == Set([
        [false false; false false],
        [false true; false false],
        [false false; true false],
        [false true; true false],
        [false false; false true],
        [false true; false true],
        [false false; true true],
        [false true; true true]])

    @test length(o.membranes) == 1
    @test length(o.solvers) == 2
    @test length(o.diags) == 2
    @test isempty(o.solvers[2])

    assignment = only(o.solvers[1].stack)
    @test isempty(assignment.clauses)
    @test isempty(assignment.free)
    @test only(assignment.affirmed) == eg.Atom(1, 1)

    eg.step!(o)
    @test isequal(eg.eugrid(o), [true missing; missing missing])

    eugrids = map(1:8) do _
        eg.step!(o)
        eg.eugrid(o)
    end
    @test Set(eugrids) == Set([
        [true false; false false],
        [true true; false false],
        [true false; true false],
        [true true; true false],
        [true false; false true],
        [true true; false true],
        [true false; true true],
        [true true; true true]])

    @test length(o.solvers) == 2
    @test all(isempty, o.solvers)

    eg.step!(o)
    @test isempty(o.solvers)
    @test isempty(o.diags)


    #@test isequal(eg.eugrid(o), [false false; false false])

    #eg.step!(o)
    #@test isequal(eg.eugrid(o), [false true; false false])


end

#=
@testset "layer" begin
    t = eg.Triple(3, 4)

    r = eg.Atom.([(1, 1), (2, 1), (1, 2)])
    l = eg.Layer(r)
    @test isempty(l.diags)
    @test sort(collect(keys(l.ds))) == r
    @test solutions(l.solver) == map(x->eg.Atom.(x), [
        [], [(1, 1)], [(1, 1), (2, 1)], [(1, 1), (2, 1), (1, 2)],
        [(1, 1), (1, 2)], [(2, 1)], [(2, 1), (1, 2)], [(1, 2)]])

    r = eg.Atom.([(2, 1), (1, 2), (2, 2)])
    l = eg.wrap!(eg.Layer([eg.Atom(1, 1)]), r, [t])
    @test isempty(l.diags)
    @test sort(collect(keys(l.ds))) == r
end
=#

#=

    @test isequal(eg.eugrid(eg.wrap!(ds1, m, membrane)), [false false; false missing])
    @test isequal(eg.eugrid(eg.wrap!(ds2, m, membrane)), [true false; false missing])
    @test isequal(eg.eugrid(eg.wrap!(ds3, m, membrane)), [false false; true missing])
    @test isequal(eg.eugrid(eg.wrap!(ds4, m, membrane)), [false true; false missing])
end
=#
#=
    r = eg.Region(t, CartesianIndex(1, 2), 3)

    @test isnothing(eg.Constraint(r, eg.DistanceMatrix(fill(1, 2, 3)), false).clause)
    @test isnothing(eg.Constraint(r, eg.DistanceMatrix(fill(1, 2, 3)), true).clause)

    @test eg.Constraint(r, eg.DistanceMatrix(fill(2, 2, 3)), false).clause == [2]
    @test eg.Constraint(r, eg.DistanceMatrix(fill(2, 2, 3)), true).clause == []

    @test eg.Constraint(r, eg.DistanceMatrix(fill(3, 2, 3)), false).clause == []
    @test eg.Constraint(r, eg.DistanceMatrix(fill(3, 2, 3)), true).clause == []
end

@testset "update_clause" begin
    t = eg.Triple(3, 4)
    r = eg.Region(t, CartesianIndex(1, 2), 3)

    c = eg.Constraint(r, eg.DistanceMatrix(fill(2, 2, 3)), false)
    eg.update_clause!(c, eg.DistanceMatrix(fill(1, 3, 3)), false)
    @test isnothing(c.clause)

    c = eg.Constraint(r, eg.DistanceMatrix(fill(2, 2, 3)), false)
    eg.update_clause!(c, eg.DistanceMatrix(fill(2, 3, 3)), false)
    @test c.clause == [2, 3]

    c = eg.Constraint(r, eg.DistanceMatrix(fill(2, 2, 3)), false)
    eg.update_clause!(c, eg.DistanceMatrix(fill(2, 3, 3)), true)
    @test c.clause == [2]
end

empty_ds(n) = [eg.DistanceMatrix([i+j for i in k-1:-1:0, j in n-1:-1:0]) for k in 1:n]

@testset "empty_ds" begin
    @test empty_ds(3) == [[2 1 0], [3 2 1; 2 1 0], [4 3 2; 3 2 1; 2 1 0]]
end

@testset "mono_cnf" begin
    t = eg.Triple(3, 4)

    @test eg.pose([t], empty_ds(1)).clauses == []
    @test eg.pose([t], empty_ds(2)).clauses == []
    @test eg.pose([t], empty_ds(3)).clauses == [Set([1, 2]), Set([1, 2, 3])]
    @test eg.pose([t], empty_ds(4)) == nothing

    cnf = eg.MonoCNF([Set([1, 2]), Set([1, 2, 3]), Set([1, 4])], Set(1:4))
    @test eg.literal_counts(cnf) == Dict(1=>3, 2=>2, 3=>1, 4=>1)

    cnf = eg.MonoCNF([Set([1, 2]), Set([1, 2, 3]), Set([1, 4])], Set(1:4))
    @test eg.simplify!(cnf).clauses == [Set([1, 2]), Set([1, 4])]
    return
    cnf = eg.MonoCNF([Set([1, 2]), Set([1, 2, 3]), Set([1, 4])], Set(1:4))
    @test eg.solve!(cnf).affirmed == [1]

    cnf = eg.MonoCNF([Set([1, 2]), Set([3])], Set(1:3))
    @test eg.solve!(cnf).affirmed == [2, 3]

    @test eg.solution([t], empty_ds(3)) == [false, true, false]
end
=#
=#
=#

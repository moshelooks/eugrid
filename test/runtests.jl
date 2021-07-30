using Eugrid
using Test

@testset "triple" begin
    t = Eugrid.Triple(3, 4)
    @test t.c == 5
    @test_throws InexactError Eugrid.Triple(1, 1)

    @test collect(Eugrid.RegionIndices(t, 1, 2)) == [CartesianIndex(1, 1)]

    @test collect(Eugrid.RegionIndices(t, 1, 3)) == CartesianIndex.([1, 2, 1], [1, 1, 2])
    @test collect(Eugrid.RegionIndices(t, 2, 3)) == [CartesianIndex(2, 2)]

    @test collect(Eugrid.RegionIndices(t, 1, 4)) == CartesianIndex.(
        [2, 3, 1, 1, 1], [1, 1, 1, 2, 3])
    @test collect(Eugrid.RegionIndices(t, 2, 4)) == CartesianIndex.(
        [2, 3, 2], [2, 2, 3])
    @test collect(Eugrid.RegionIndices(t, 3, 4)) == [CartesianIndex(3, 3)]

    @test collect(Eugrid.RegionIndices(t, 1, 5)) == CartesianIndex.(
        [3, 4, 1, 1, 1], [1, 1, 2, 3, 4])
    @test collect(Eugrid.RegionIndices(t, 2, 5)) == CartesianIndex.(
        [3, 4, 2, 2, 2], [2, 2, 2, 3, 4])
    @test collect(Eugrid.RegionIndices(t, 3, 5)) == CartesianIndex.(
        [3, 4, 3], [3, 3, 4])
    @test collect(Eugrid.RegionIndices(t, 4, 5)) == [CartesianIndex(4, 4)]

    t = Eugrid.Triple(6, 8)
    @test t.c == 10

    @test collect(Eugrid.RegionIndices(t, 1, 8)) == CartesianIndex.(
        [3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7])

    @test collect(Eugrid.RegionIndices(t, 2, 8)) == CartesianIndex.(
        [3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 2],
        [2, 2, 2, 2, 2, 1, 2, 3, 4, 5, 6, 7])
end

@testset "region" begin
    t = Eugrid.Triple(3, 4)

    r = Eugrid.Region(t, CartesianIndex(1, 1), 1)
    @test r.v == r.h == 1:1
    @test Eugrid.lastk(r) == 1
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 5

    r = Eugrid.Region(t, CartesianIndex(1, 1), 2)
    @test r.v == r.h == 1:2
    @test Eugrid.lastk(r) == 2
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 5

    r = Eugrid.Region(t, CartesianIndex(1, 1), 3)
    @test r.v == r.h == 1:3
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 5

    r = Eugrid.Region(t, CartesianIndex(1, 1), 4)
    @test r.v == 1:3
    @test r.h.start == 1
    @test length(r.h) == 0
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 5

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(1, 1), 5)

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(2, 1), 1)

    r = Eugrid.Region(t, CartesianIndex(2, 1), 2)
    @test r.v == 2:2
    @test r.h == 1:2
    @test Eugrid.lastk(r) == 2
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 5

    r = Eugrid.Region(t, CartesianIndex(2, 1), 3)
    @test r.v == 2:3
    @test r.h == 1:3
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 5

    r = Eugrid.Region(t, CartesianIndex(2, 1), 4)
    @test r.v == 2:4
    @test r.h == 1:4
    @test Eugrid.lastk(r) == 4
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 5

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(2, 1), 5)

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(1, 2), 1)

    r = Eugrid.Region(t, CartesianIndex(1, 2), 2)
    @test r.v == 1:2
    @test r.h == 2:2
    @test Eugrid.lastk(r) == 2
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 6

    r = Eugrid.Region(t, CartesianIndex(1, 2), 3)
    @test r.v == 1:3
    @test r.h == 2:3
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 6

    r = Eugrid.Region(t, CartesianIndex(1, 2), 4)
    @test r.v == 1:3
    @test r.h.start == 2
    @test length(r.h) == 0
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 6

    r = Eugrid.Region(t, CartesianIndex(1, 2), 5)
    @test r.v == 1:3
    @test r.h.start == 2
    @test length(r.h) == 0
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 4
    @test Eugrid.lastcol(r) == 6

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(1, 2), 6)

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(2, 2), 1)

    r = Eugrid.Region(t, CartesianIndex(2, 2), 2)
    @test r.v == r.h == 2:2
    @test Eugrid.lastk(r) == 2
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 6

    r = Eugrid.Region(t, CartesianIndex(2, 2), 3)
    @test r.v == r.h == 2:3
    @test Eugrid.lastk(r) == 3
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 6

    r = Eugrid.Region(t, CartesianIndex(2, 2), 4)
    @test r.v == r.h == 2:4
    @test Eugrid.lastk(r) == 4
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 6

    r = Eugrid.Region(t, CartesianIndex(2, 2), 5)
    @test r.v == 2:4
    @test r.h.start == 2
    @test length(r.h) == 0
    @test Eugrid.lastk(r) == 4
    @test Eugrid.lastrow(r) == 5
    @test Eugrid.lastcol(r) == 6

    @test_throws AssertionError Eugrid.Region(t, CartesianIndex(2, 2), 6)
end

@testset "regions" begin
    t = Eugrid.Triple(3, 4)
    @test collect(Eugrid.regions(t, 3, 4)) == Eugrid.Region.(
        [t], Eugrid.RegionIndices(t, 3, 4), 4)

end

@testset "negation" begin
    t = Eugrid.Triple(3, 4)

    r = Eugrid.Region(t, CartesianIndex(1, 1), 3)
    @test Eugrid.negation_distance(r, 3, 3) == 2

    r = Eugrid.Region(t, CartesianIndex(1, 1), 4)
    @test Eugrid.negation_distance(r, 3, 4) == 3

    r = Eugrid.Region(t, CartesianIndex(2, 1), 4)
    @test Eugrid.negation_distance(r, 3, 4) == 2
    @test Eugrid.negation_distance(r, 4, 4) == 3

    d = Eugrid.DistanceMatrix(fill(2, 3, 4))
    @test Eugrid.negates(r, d)
    d[1, 2] = d[2, 1] = 3
    @test !Eugrid.negates(r, d)

    r = Eugrid.Region(t, CartesianIndex(1, 2), 3)
    @test Eugrid.negation_distance(r, 1, 3) == -1
    d = Eugrid.DistanceMatrix(fill(-1, 1, 3))
    @test Eugrid.negates(r, d)
    d[1, 2] = 2
    @test !Eugrid.negates(r, d)
end

@testset "affirmation" begin
    t = Eugrid.Triple(3, 4)

    r = Eugrid.Region(t, CartesianIndex(1, 1), 2)
    @test Eugrid.affirmation_distance(r, 2, 1) == 1
    @test Eugrid.affirmation_distance(r, 2, 2) == 2
    @test Eugrid.affirmation_distance(r, 1, 2) == 2

    d = Eugrid.DistanceMatrix([2 9; 9 9])
    @test Eugrid.affirmation_bound(r, d) == 0

    r = Eugrid.Region(t, CartesianIndex(2, 1), 3)
    @test Eugrid.affirmation_distance(r, 3, 1) == 1
    @test Eugrid.affirmation_distance(r, 3, 2) == 2
    @test Eugrid.affirmation_distance(r, 3, 3) == 3
    @test Eugrid.affirmation_distance(r, 2, 3) == 2

    r = Eugrid.Region(t, CartesianIndex(2, 1), 4)
    @test Eugrid.affirmation_distance(r, 4, 2) == 2
    @test Eugrid.affirmation_distance(r, 4, 3) == 3
    @test Eugrid.affirmation_distance(r, 4, 4) == 4
    @test Eugrid.affirmation_distance(r, 1, 4) == 1
    @test Eugrid.affirmation_distance(r, 2, 4) == 2
    @test Eugrid.affirmation_distance(r, 3, 4) == 3

    r = Eugrid.Region(t, CartesianIndex(1, 2), 3)
    @test Eugrid.affirmation_distance(r, 3, 2) == 1
    @test Eugrid.affirmation_distance(r, 3, 3) == 2
    @test Eugrid.affirmation_distance(r, 2, 3) == 2
    @test Eugrid.affirmation_distance(r, 1, 3) == 2

    d = Eugrid.DistanceMatrix(fill(-99, 1, 3))
    d[1, 2] = 5
    @test Eugrid.affirmation_bound(r, d) == 3

    d = Eugrid.DistanceMatrix(fill(-99, 2, 3))
    d[1, 2] = 1
    @test Eugrid.affirmation_bound(r, d) == -1

    d[1, 2] = 2
    @test Eugrid.affirmation_bound(r, d) == -100

    d[1, 2] = 9
    d[2, 1] = 99
    @test Eugrid.affirmation_bound(r, d) == 7
end

@testset "constraint" begin
    t = Eugrid.Triple(3, 4)
    r = Eugrid.Region(t, CartesianIndex(1, 2), 3)

    @test isnothing(Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(1, 2, 3)), false).clause)
    @test isnothing(Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(1, 2, 3)), true).clause)

    @test Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(2, 2, 3)), false).clause == [2]
    @test Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(2, 2, 3)), true).clause == []

    @test Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(3, 2, 3)), false).clause == []
    @test Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(3, 2, 3)), true).clause == []
end

@testset "update_clause" begin
    t = Eugrid.Triple(3, 4)
    r = Eugrid.Region(t, CartesianIndex(1, 2), 3)

    c = Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(2, 2, 3)), false)
    Eugrid.update_clause!(c, Eugrid.DistanceMatrix(fill(1, 3, 3)), false)
    @test isnothing(c.clause)

    c = Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(2, 2, 3)), false)
    Eugrid.update_clause!(c, Eugrid.DistanceMatrix(fill(2, 3, 3)), false)
    @test c.clause == [2, 3]

    c = Eugrid.Constraint(r, Eugrid.DistanceMatrix(fill(2, 2, 3)), false)
    Eugrid.update_clause!(c, Eugrid.DistanceMatrix(fill(2, 3, 3)), true)
    @test c.clause == [2]
end

empty_ds(n) = [Eugrid.DistanceMatrix([i+j for i in k-1:-1:0, j in n-1:-1:0]) for k in 1:n]

@testset "empty_ds" begin
    @test empty_ds(3) == [[2 1 0], [3 2 1; 2 1 0], [4 3 2; 3 2 1; 2 1 0]]
end

@testset "mono_cnf" begin
    t = Eugrid.Triple(3, 4)

    @test Eugrid.MonoCNF([t], empty_ds(1)).clauses == []
    @test Eugrid.MonoCNF([t], empty_ds(2)).clauses == []
    @test Eugrid.MonoCNF([t], empty_ds(3)).clauses == [Set([1, 2]), Set([1, 2, 3])]
    @test_throws AssertionError Eugrid.MonoCNF([t], empty_ds(4))

    cnf = Eugrid.MonoCNF([Set([1, 2]), Set([1, 2, 3]), Set([1, 4])])
    @test Eugrid.literal_counts(cnf) == Dict(1=>3, 2=>2, 3=>1, 4=>1)

    cnf = Eugrid.MonoCNF([Set([1, 2]), Set([1, 2, 3]), Set([1, 4])])
    @test Eugrid.simplify!(cnf).clauses == [Set([1, 2]), Set([1, 4])]

    cnf = Eugrid.MonoCNF([Set([1, 2]), Set([1, 2, 3]), Set([1, 4])])
    @test Eugrid.solve!(cnf).clauses == [Set([1])]

    cnf = Eugrid.MonoCNF([Set([1, 2]), Set([3])])
    @test Eugrid.solve!(cnf).clauses == [Set([2]), Set([3])]

    @test Eugrid.solution([t], empty_ds(3)) == [false, true, false]
end

#=








julia>
Eugrid.Region(Eugrid.Triple(3, 4, 5), 2:4, 2:1, false, Int64[])

@testset "triple_indices" begin

    @test collect(Eugrid.TripleIndices(1, 1, 3, 4)) == CartesianIndex.([], [])

    @test collect(Eugrid.TripleIndices(2, 1, 3, 4)) == CartesianIndex.([], [])

    @test collect(Eugrid.TripleIndices(2, 2, 3, 4)) == CartesianIndex.([1], [1])

    @test collect(Eugrid.TripleIndices(3, 1, 3, 4)) == CartesianIndex.([], [])

    @test collect(Eugrid.TripleIndices(3, 2, 3, 4)) == CartesianIndex.([1, 2], [1, 1])

    @test collect(Eugrid.TripleIndices(3, 3, 3, 4)) == CartesianIndex.(
        [1, 2, 1, 2],
        [1, 1, 2, 2])

    @test collect(Eugrid.TripleIndices(4, 1, 3, 4)) == CartesianIndex.(
        [],
        [])

    @test collect(Eugrid.TripleIndices(4, 2, 3, 4)) == CartesianIndex.(
        [1, 2, 3],
        [1, 1, 1])

    @test collect(Eugrid.TripleIndices(4, 3, 3, 4)) == CartesianIndex.(
        [1, 2, 3, 1, 2, 3],
        [1, 1, 1, 2, 2, 2])

    @test collect(Eugrid.TripleIndices(4, 4, 3, 4)) == CartesianIndex.(
        [1, 2, 3, 1, 2, 3, 1, 2],
        [1, 1, 1, 2, 2, 2, 3, 3])

    @test collect(Eugrid.TripleIndices(5, 1, 3, 4)) == CartesianIndex.(
        [],
        [])

    @test collect(Eugrid.TripleIndices(5, 2, 3, 4)) == CartesianIndex.(
        [1, 2, 3],
        [1, 1, 1])

    @test collect(Eugrid.TripleIndices(5, 3, 3, 4)) == CartesianIndex.(
        [1, 2, 3, 1, 2, 3],
        [1, 1, 1, 2, 2, 2])

    @test collect(Eugrid.TripleIndices(5, 4, 3, 4)) == CartesianIndex.(
        [1, 2, 3, 1, 2, 3, 1, 2],
        [1, 1, 1, 2, 2, 2, 3, 3])

    @test collect(Eugrid.TripleIndices(5, 5, 3, 4)) == CartesianIndex.(
        [1, 2, 3, 1, 2, 3, 1, 2],
        [1, 1, 1, 2, 2, 2, 3, 3])

    @test collect(Eugrid.TripleIndices(5, 5, 3, 5)) == CartesianIndex.(
        [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 1, 2],
        [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4])

    @test collect(Eugrid.TripleIndices(5, 4, 3, 5)) == CartesianIndex.(
        [1, 2, 3, 4, 1, 2, 3, 4, 1, 2],
        [1, 1, 1, 1, 2, 2, 2, 2, 3, 3])

    @test collect(Eugrid.TripleIndices(5, 4, 2, 6)) == CartesianIndex.(
        [1, 2, 3, 4, 1, 1],
        [1, 1, 1, 1, 2, 3])

end


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
=#

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
    @test_throws Eugrid.UnsatisfiableException Eugrid.MonoCNF([t], empty_ds(4))

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

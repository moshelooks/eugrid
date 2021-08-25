using Eugrid
using Test

const eg = Eugrid

@testset "distance_matrix" begin
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
#=
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

@testset "constrain!" begin
    t = eg.Triple(3, 4)

    cs = eg.Constraints()
    d = eg.diag(eg.DistanceMatrix(eg.Atom(1, 1)))
    eg.constrain!(cs, eg.Box([t]), d)
    @test cs.domain == Set([eg.Atom(2, 2)])
    @test cs.region_clauses == Dict(
        eg.Region(t, eg.Atom(1, 1))=>eg._free,
        eg.Region(t, eg.Atom(2, 1))=>eg._free,
        eg.Region(t, eg.Atom(1, 2))=>[eg.Atom(2, 2)],
        eg.Region(t, eg.Atom(2, 2))=>eg._free)

    cs = eg.Constraints()
    d = eg.diag(d)
    eg.constrain!(cs, eg.Box([t]), d)
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
    eg.constrain!(cs, eg.Box([t]), d)
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
    ds = filter(square(Set([eg.Atom(1, 1), eg.Atom(2, 2)]), 3)) do (a, _)
        maximum(a.I) == 3
    end
    cs = eg.constraints(eg.Box([t]), ds)
    @test cs.domain == Set(eg.Atom.([(3, 1), (3, 2), (1, 3), (2, 3)]))
    expected = Dict(eg.Region(t, a)=>eg._free for a in eg.Atoms(3))
    @test cs.region_clauses == expected

    ds = filter(square(Set([eg.Atom(2, 2)]), 3)) do (a, _)
        maximum(a.I) == 3
    end
    cs = eg.constraints(eg.Box([t]), ds)
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
    @test o.box.ts == eg.Triple.([4, 3], [3, 4])
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
=#



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

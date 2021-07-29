struct Triple
    a::Int
    b::Int
    c::Int
    Triple(a::Int, b::Int) = new(a, b, sqrt(a^2 + b^2))
end

vrange(t::Triple, m::Int, k::Int) = (k + t.b > m ? k : max(k, m-t.a+1)):m-1
hrange(t::Triple, m::Int, k::Int) = max(k+1, m-t.b+1):m-1

RegionIndices(t::Triple, m::Int, k::Int) =
    Iterators.flatten((CartesianIndices((vrange(t, m, k), k:k)),
                       CartesianIndices((k:k, hrange(t, m, k)))))

struct Region
    t::Triple
    v::UnitRange{Int}
    h::UnitRange{Int}

    function Region(t::Triple, u::CartesianIndex{2}, m::Int)
        vstart, hstart = u.I
        @assert max(vstart, hstart) <= m
        vstop, hstop = min.((vstart, hstart) .+ (t.a, t.b) .- 1, m)
        @assert max(vstop, hstop) >= m
        if vstop < m
            hstop = hstart-1
        elseif hstop < m
            vstop = vstart-1
        end
        new(t, vstart:vstop, hstart:hstop)
    end
end

vorigin(r::Region)::CartesianIndex{2} = CartesianIndex(r.h.start, r.v.start)
horigin(r::Region)::CartesianIndex{2} = CartesianIndex(r.v.start, r.h.start)

lastrow(r::Region)::Int = r.v.start + r.t.a
lastcol(r::Region)::Int = r.h.start + r.t.b

lastk(r::Region)::Int =
    isempty(r.v) ? r.h.stop : isempty(r.h) ? r.v.stop : max(r.h.stop, r.v.stop)

regions(t::Triple, m::Int, k::Int) = (Region(t, u, m) for u in RegionIndices(t, m, k))

negation_distance(r::Region, m::Int, k::Int)::Int = r.t.c + m + k - lastrow(r) - lastcol(r)

function negates(r::Region, d::Matrix{Int})::Bool
    m, k = size(d)
    target = negation_distance(r, m, k)
    k in r.h && d[horigin(r)] == target && return true
    k in r.v && d[vorigin(r)] == target && return true
    false
end

affirmation_distance(r::Region, m::Int, k::Int)::Int =
    r.t.c - max(lastrow(r) - m, lastcol(r) - k)

function affirmation_bound(r::Region, d::Matrix{Int})::Int
    m, k = size(d)
    if k in r.h
        bound = d[horigin(r)] - affirmation_distance(r, m, k)
        if bound >= 0 && k in r.v
            bound = min(bound, d[vorigin(r)] - affirmation_distance(r, k, m))
        end
    else
        @assert k in r.v
        bound = d[vorigin(r)] - affirmation_distance(r, k, m)
    end
    bound
end

mutable struct Constraint
    region::Region
    clause::Union{Vector{Int}, Nothing}
end

function Constraint(r::Region, d::Matrix{Int}, negated::Bool)
    bound = affirmation_bound(r, d)
    bound < 0 && return Constraint(r, nothing)
    c = Vector{Int}()
    bound == 0 && !negated && push!(c, size(d)[2])
    Constraint(r, c)
end

function update_clause!(c::Constraint, d::Matrix{Int}, negated::Nothing)::Bool
    !negates(c.region, d) && return false
    c.clause = nothing
    true
end

function update_clause!(c::Constraint, d::Matrix{Int}, negated::Bool)::Nothing
    isnothing(c.clause) && return
    bound = affirmation_bound(c.region, d)
    if  bound < 0
        c.clause = nothing
    elseif bound == 0 && !negated
        push!(c.clause, size(d)[2])
    end
    nothing
end

constraints(ts::Vector{Triple}, d::Matrix{Int}, negated::Bool) =
    (Constraint(r, d, negated) for t in ts for r in regions(t, size(d)...))

struct MonoCNF
    clauses::Vector{Set{Int}}
end

function MonoCNF(ts::Vector{Triple}, ds::Vector{Matrix{Int}})::MonoCNF
    clauses = Vector{Vector{Int}}()
    active_constraints = Vector{Constraint}()
    for (k, d) in enumerate(ds)
        negated = mapreduce(|, active_constraints, init=false) do c
            update_clause!(c, d, nothing)
        end
        filter!(active_constraints) do c
            update_clause!(c, d, negated)
            k != lastk(c.region) && return true
            if !isnothing(c.clause)
                @assert !isempty(c.clause)
                push!(clauses, c.clause)
            end
            false
        end
        append!(active_constraints, constraints(ts, d, negated))
    end
    MonoCNF(Set{Int}.(sort!(clauses)))
end

literal_counts(cnf::MonoCNF)::Dict{Int, Int} =
    mapreduce(c->Dict(c.=>1), mergewith(+), cnf.clauses)

function simplify!(cnf::MonoCNF)::MonoCNF
    for subset in sort(cnf.clauses, by=length)
        filter!(cnf.clauses) do c
            !(length(subset) < length(c) && issubset(subset, c))
        end
    end
    cnf
end

function solve!(cnf::MonoCNF)::MonoCNF
    isempty(cnf.clauses) && return cnf
    unique!(cnf.clauses)
    while true
        simplify!(cnf)
        n, l = maximum(reverse, literal_counts(cnf))
        n == 1 && break
        push!(cnf.clauses, Set((l,)))
    end
    for c in cnf.clauses
        filter!(>=(maximum(c)), c)
    end
    sort!(cnf.clauses, by=only)
    cnf
end

function solution(ts::Vector{Triple}, ds::Vector{Matrix{Int}})::BitVector
    diags = falses(length(ds))
    diags[map(only, solve!(MonoCNF(ts, ds)).clauses)] .= true
    diags
end


#=
update_term(::Free









1. satisfied (1 or more blocks)
2. free

    clause::


function update_clause!(r::Region, k::Int, bound::Int)::Bool
    bound > r.t.c && return false
    if bound == r.t.c
        r.clause = nothing
    else
        push!(r.clause, k)
    end
    true
end

function update_clause!(r::Region, d::Matrix{Int})::Nothing
    isnothing(r.clause) && return
    m, k = size(d) .+ 1
    bound = r.t
    if k in r.v
        vbound = d[v.start, h.start] + max(lastrow(r) - k, lastcol(r) - m)
        if vbound == r.t.c
            r.clause = nothing
            return
        end
    end
    if k in r.h
        hbound = d[h.start, v.start] + max(lastrow(r) - m, lastcol(r) - k)
        if hbound == r.t.c
            r.clause = nothing
            return
        end
    end
    ((k in r.v && vbound < r.t.c) || (k in r.h && hbound < r.t.c)) && push!(r.clause, k)
end


CNF = Vector{Set{Int}}

mut

function build_cnf(ts::Vector{Triple}, ds::Vector{Matrix{Int}})::CNF
    cnf = CNF()
    active_regions = Vector{Region}()
    for (k, d) in enumerate(ds)
        blocked = false
        filter!(active_regions) do r
            blocked = blocked || blocks(r, d)
            k == lastk(r) && return true
            if !isnothing(r.clause)
                @assert !isempty(r.clause)
                push!(cnf, Set(r.clause))
            end
            false
        end
        append!(active_regions, regions(ts, k))
        blocked && continue
        for r in active_regions
            update_clause!(r, k)
        end
    end
    cnf
end



    target = r.t.c + sum(size(d)) + 2
    isbottom(r, d) && d[r.tl] == target  && return true
    isright(r, d) && d[transpose(r.tl)] == target && return true
    return false
end
    target = r.c + sum(mk(d)) - sum(r.br.I)



struct ClauseBuilder
    triples::Vector{Triple}
    atom_distances::Vector{Matrix{Int}}
end


TripleIndices(m::Int, k::Int, a::Int, b::Int) =
    Iterators.flatten((CartesianIndex(1,1):CartesianIndex(min(b,m)-1,min(a,k)-1),
                       CartesianIndex(1,a):CartesianIndex(min(a,m)-1,min(b,k)-1)))



function regions(cb::ClauseBuilder, k::Int)
    m = length(cb.ds) + 1
    k == m && return []
    for t in cb.triples
        for ul in CartesianIndex()
        end
    end
end

islast(r::Region, k::Int, m::Int)::Int = k == m || k == minimum(r.br.I)

function isbottom(r::Region, d::DMatrix)::Bool
    m, k = mk(d)
    lastrow(r) > m &&





    r.ul[2] < k < r.lr[2] && d[r.ul] + sum(r.lr.I) - m - k == r.t.c && return true
    r.

    k = CartesianIndex(size(d) .+ 1)
    k

    d[r.ul



            r.satisfied && continue
            bound = satisfaction_bound(r, k)
            if bound < 0
                r.satisfied = true
            elseif bound == 0
                push!(r.clause, k)
            end



    end

        for ul in CartesianIndex(m-1,k):CartesianIndex(m-t.a+1,k)


function distance(ws::Workspace, r::Region)
    i[1] == r.ul[1] && return
end

function distance(ws::Workspace, u::CartesianIndex{2}, v::CartesianIndex{2})
    v[1] == u[1] && return v[2]-u[2]
    v[2] == u[2] && return v[1]-u[1]
    v[1] == m && return ws.d[u][v[2]-


function clause(ws::Workspace, r::Region)::Union[Nothing, Set{Int}]
    c = Vector{Int}()
    for i in ribbon(ws, r)
        blocked(ws, i) && continue
        d = distance(ws, r.ul, i) + maximum((r.lr - i).I)
        d < r.t.c && return nothing
        d == r.t.c && push!(c, minimum(i.I))
    end
    Set(c)


    if r.lr[1] > m
        va
    for (d, i) in ribbon(ws, r)
        d +=



function ribbon(r::Region)
    iters = []
end

Base.iterate(r::Ribbon, u::CartesianIndex{2}=r.first) =
    (r.



struct Cell
    u::CartesianIndex
    d::Vector{Int}
end

function grow!(c::Cell, diags::BitVector)::Nothing



struct Triple
    a::Int
    b::Int
    c::Int
end

struct Region
    ul::CartesianIndex
    t::Triple
end

function available(r::Region, d::Matrix{Int})::Bool
    i = CartesianIndex(size(d) .+ 1)
end



regions(triples::Vector{Triple}, m::Int)::Vector{Region}

Clause = Vector{Int}

mutable struct Workspace
    triples::Vector{Triple}
    available::BitVector
    regions::Vector{Region}
    cnf::Vector{Clause}
    counts::Dict{Int, Int}
end

Workspace(triples::Vector{Triple}, m::Int) =
    Workspace(triples, trues(m), regions(triples, m), [], [])

function foo(ws::Workspace)
    for k in findall(ws.available)
        for r in ws.regions






struct PinnedTriple
    p::CartesianIndex{2}
    t::Triple
end

Base.CartesianIndices(pt::PinnedTriple) = p:p+CartesianIndex(t)

TripleIndices(

struct Literal
    d::Matrix{Int}
end

function grow(







isblocked(ts::Vector{Triple}, d::Matrix{Int}) =
    any(i->sum(i.I))


struct Border
    ts::Vector{Triple}
    ds::Vector{Matrix{Int}}
    m::Int
    blocked::Set{Int}
end

Border(ts::Vector{Triple}, ds::Vector{Matrix{Int}}) =
    Border(ts, ds, length(ds) + 1, filter(i->isblocked(ts, ds[i]), eachindex(ds)))


function fringe(e::Edge, t::Triple, u::CartesianIndex{2})::Vector{KDAB}
    i
    a1 = border.m - u[1]
    if a1 < a
        a2 = a - a1
        push!(kdabs, (u[2], a1, a2, b))
        for (i, k) in enumerate(u[2]+1:min(u[2]+b, m))
            f(k, d[k][a1, i], a2, b - i)
        end
    end

    b1 = m - u[2]
    if b1 < b
        b2 = b - b1
        f(u[1], b1, a, b2)
        for (i, k) in enumerate(u[1]+1:min(u[1]+a, m - 1))
            f(k, d[k][b1, i], a - i, b2)
        end
    end
end



function findclause(f::Fringe, t::Triple, u::CartesianIndex{2})
    for (k, d, v) in fringe(f, t, u)


clauses(f::Fringe) =
    filter(!isnothing, (findclause(f, t, u) for t in f.ts for u in DistanceIndices(f.m, t)))
        clause = findclause(t
        !isnothing(clause) && push!(clauses, clause)
    end
    Constraints(ts, ds, blocked, clauses)


    union!(f.blocked, filter(i->o


TripleIndices(m::Int, k::Int, a::Int, b::Int) =
    Iterators.flatten((CartesianIndex(1,1):CartesianIndex(min(b,m)-1,min(a,k)-1),
                       CartesianIndex(1,a):CartesianIndex(min(a,m)-1,min(b,k)-1)))

isblocked(d::Matrix{Int}, a::Int, b::Int, c::Int)::Bool =
    any(i->sum(i.I) - d[i] == c, TripleIndices((size(d) .+ 1)..., a, b))

isblocked(d::Matrix{Int}, ts::Vector{Pythagorean.Triple})::Bool =
    any(t->isblocked(d, t...), ts)

isblocked(ds::Vector{Matrix{Int}}, ts::Vector{Pythagorean.Triple})::Set{Int} =
    Set{Int}(i for i in eachindex(ds) if isblocked(ds[i], ts))

DistanceIndices(m::Int, a::Int, b::Int) =
    Iterators.flatten((CartesianIndex(max(1, m-a+1),1):CartesianIndex(m-1, m-b),
                       CartesianIndex(1, max(1, m-b+1)):CartesianIndex(m-1, m-1)))


#function fringe(::Int, u::CartesianIndex, a::Int, b::Int)
#    iters = []

function findclause(t::Pythagorean.Triple, ds::Vector{Matrix{Int}}, u::CartesianIndex)
    mapfringe()
    for
    end
end



mutable struct Square
    ts::Vector{Pythagorean.Triple}
    ds::Vector{Matrix{Int}}
    blocked::Set{Int}
    clauses::Vector{Set{Int}}
end

function Constraints(ts::Vector{Pythagorean.Triple}, ds::Vector{Matrix{Int}})::Constraints
    m = length(ds) + 1
    blocked = Set{Int}(i for i in 1:m if isblocked(ts, ds))
    clauses = Vector{Set{Int}}()
    for (a, b, c) in ts, u in DistanceIndices(m, a, b)
        clause = findclause(ds, a, b, c)
        !isnothing(clause) && push!(clauses, clause)
    end
    Constraints(ts, ds, blocked, clauses)
end

function grow(n::Int)::AbstractMatrix{Bool}
    diags = BitMatrix(n)
    constraints = Constraints(Pythagorean.triples(n))
    for i in 1:n
        diags[1:i, i] .= satisfy!(constraints)
    state = State(prev)



function DistanceIndices(f, m::Int, u::CartesianIndex{2}, a::Int, b::Int)
    a1 = m - u[1]
    if a1 < a
        a2 = a - a1
        f(u[2], a1, a2, b)
        for (i, k) in enumerate(u[2]+1:min(u[2]+b, m))
            f(k, d[k][a1, i], a2, b - i)
        end
    end

    b1 = m - u[2]
    if b1 < b
        b2 = b - b1
        f(u[1], b1, a, b2)
        for (i, k) in enumerate(u[1]+1:min(u[1]+a, m - 1))
            f(k, d[k][b1, i], a - i, b2)
        end
    end
end

Clause = Set{Int}

function clause(ds::Vector{Matrix{Int}}, t::Triple,
                blocked::Set{Int})::Union[Clause, Nothing]
    c = Clause()
    DistanceIndices(m, u, a, b) do (k, d, a2, b2)
        k in blocked && continue
        dmin = d + max(a2, b2)
        dmin < c && return nothing
        dmin == t[3] & push!(c, k)
    end
    c
end

function clauses(ds::Vector{Matrix{Int}}, ts::Vector{Triple},
                 blocked::Set{Int})::Vector{Clause}
    for (a, b, c) in ts, u in DistanceIndices(
        c = clause(ds, a, b, c, blocked)
        isnothing(c) && continue
        push!(c, clauses)




    for i in 1:m-1
        for (a, b, c) in ts
            for (k, j) in DistanceIndices(i, a, b)
                k in blocked && continue
                dmin = sum(j.I) - d[k][j]

            for i in ClauseIndices(k, a, b)
                for j in DistanceIndices(k,

            for j in 1:min(a,m


    end

         && return true




    d = s.d[i]
    for (a, b, c) in s.triples
        bound = a + b - c
        for y in 1:min(a,i)-1, x in 1:min(b,m)-1
            d[x
        end
        for y in a:min(b,i)-1, x in 1:min(a,m)-1



        for y in 1:min(b,m)-1, x in 1:min(a,y <= i ? m : i)-1
            x + y - s.d[i][x, y] == bound && return true
        end
    end
    false
end



        for y in 1:min(b,i)-1, x in 1:min(a,m)-1
            x + y - s.d[i][x, y] == bound && return true
        end
        for y in 1:min(a,i)-1, x in min(b,i):min(b,m)-1
            x + y - s.d[i][x, y]
    end
    false
end



        for y in max(1, i - b + 1):i-1
            for x in max(1, m - a + 1):m-1
                s.d[


function foo()
    for i in 1:m
        for (d, a, b, c) in foo(i)
        end
    end
end

function bar()
    for (u, a, b, c) in triples(m)
        for

function regions()
    for ul in CartesianIndices((m-1, m-1))
        for t in triples(m - ul[1], m - ul[2])



function isblocked(d::Matrix{Int}, triples::Vector{Triple})
    for i in CartesianIndices(d), (a, b, c) in triples
        (i[1] >= a || i[2] >= b) && continue
        d[i] + a + b - sum(i.I) == c && return true
    end
    return false
end

function grow!(ds::Vector{Matrix{Int}}, triples::Vector{Triple})
    blocked = Set(filter(i->isblocked(d[i], triples), 1:length(ds)))
    for (a, b, c) in triples
        for i in 1:n-





struct Square
    n::Int
    d::Array{Int, 4}
    diags::CircularArray{Bool, 2, BitMatrix}
    triples::Matrix{Vector{Tuple{CartesianIndex{2}, Int}}}
end

function distances(t, u, d)
    dy = d[u[2]]
    n, m = size(dy) .+ 1
    d = dy[1:min(b, m)-1, n]
    a2 = a - (n - u[2]) + 1
    dlower = d .+ [max(b - i + 1, a2) for i in eachindex(d)]
    dupper = d .+ [(b - i + a2) for i in eachindex(d)]
    for i in eachindex(dupper)
        if dupper[i] == t.c
            push!(blocked, u[2] + i)


function grow()
    n, m = 1 .+ size(d)
    d = [[(1 .+ d) (2 .+ view(d, :, m - 1))]; (2:m+1)']


Triple = Tuple{Int, Int, Int}

struct ROI
    t::Triple
    u::CartesianIndex{2}
end

function foo()
    for roi
        for duv in d[roi.u[2]][2:end, roi.u[1]]


function all_rois(n::Int, m::Int)
    @assert n <= m
    candidates = ((a, b, sqrt(a^2 + b^2)) for b in 1:m for a in 1:b)
    triples = [(a, b, Int(c)) for (a, b, c) in candidates if isinteger(c)]
    rois = Vector{ROI}()
    for t in triples
        a, b, _ = t
        start = CartesianIndex(max(1, n - a + 1), 1)
        stop = CartesianIndex(min(n, m - a + 1), min(n, m - b + 1))
        append!(rois, [ROI(t, i) for i in start:stop])
    end
    rois
end



function Square(n::Int)::Square
    d = repeat(CartesianIndices((k, k)) .|> x->sum(x.I), outer=(1, 1, n, n))
end

function constraints(f, s::Square, d::Distances)
    for (a, b, c) in s.triples
        for i in 1:length(d)
            for j in max(1, length(d) - a):length(d)
                f(i, d)


mse(s::Square)::Float64 =
    mean((s.d[a, b, i] - c)^2 for i in CartesianIndices(s) for (a, b, c) in s.triples[i])

function distances(diags::AbstractVector{Bool}, dp::Matrix{Int})
    n, m = 1 .+ size(dp)
    @assert length(diags) == m > 1
    d = Matrix{Int}(undef, n, m)
    if diags[1]
        d[:, 1] .= n:-1:1
    else
        d[1:n-1, 1] .= view(dp, :, 1) .+ 1
        d[n, 1] = d[n-1,1] + 1
    end
    for i in 2:m
        if diags[i]
            d[1:n-1, i] .= 1 .+ view(dp, :, i - 1)
            d[n, i] = i
        elseif i < m
            d[1:n-1, i] .= 1 .+ min.(view(dp, :, i), view(d, 1:n-1, i - 1))
            d[n, i] = 1 + d[n, i - 1]
        else
            d[:, i] = 1 .+ view(d, :, i - 1)
        end
    end
    d
end

function distances(diags::AbstractVector{Bool}, dp::Vector{Matrix{Int}})::Vector{Matrix{Int}}
    n = length(diags)
    @assert length(dp) == n - 1
    ds = Vector{Matrix{Int}}(undef, n)
    ds[n] = Matrix{Int}(undef, n, 1)
    for (i, dpi) in enumerate(dp)
        ds[i] = distances(view(diags, i:n), dpi)
        ds[n][i] = ds[i][n, n - i + 1]
    end
    ds[n][n] = diags[n] ? 1 : 2
    ds
end

function distances(s::Square)::Vector{Matrix{Int}}
    d = Vector{Matrix{Int}}()
    for i in 1:LinearAlgebra.checksquare(s.diags)
        d = distances(view(s.diags, 1:i, i), d)
    end
end

function grow(dp::Vector{Matrix{Int}})
    n = length(dp) + 1
    for (a, b, c) in triples(n)
        for (i,di) in enumerate(d)
            dwithout = d + b - i + 1
            @assert dwithout <= c (dwithout, c)
            if dwithout == c
                push!(blocked, ix)
            end
        end
    end
    cnf = [setdiff(ix+1:ix+b-1, blocked) for (ix, b, _, _) in triples]
    for c in cnf
        @assert !isempty(c)
    end
    diags = falses(n)
    singletons = unique(only(c) for c in cnf if length(c) == 1)
    filter!(c->!any(in(c), singletons), cnf)
    diags[singletons] = true
    while !isempty(cnf)
        l = maximum(domain, l->(count[l], l))
        diags[l] = true
        true!(cnf, l)
    end


    while true
        isempty(disjuncts) && return diags
        toremove = unique(only(dj) for dj in disjuncts if length(dj) == 1)
        isempty(toremove) && break
        diags[toremove] .= true
        for dj in disjuncts
            setdiff!(dj, toremove)
        end
        filter!(!isempty, disjuncts)

        sort!(disjuncts, by=length, rev=true)
        i = findlast(x->length(x)==1, disjuncts)
        isnothing(i) && break
        toremove = only.(disjuncts[i:end]
        resize!(disjuncts, i-1)
        for disjunct in disjuncts
            setdiff!(distjunct, toremove)

        if length(
        for ix in ix
            ix in blocked &&

        clause =
        m = n-b+1
        for d in view(ds, 1:m)
            for j in 1:m-b+1

                if d[n-
                for k in 1:
                dt = d[n-a+1,b]
                if dt == c - 2
                    push!(blocked,
                @assert c - 2 <= dt <=




    if diags[m]
        d[1:n-1, m]

    if k == 1
        d[:, 1] .= n:-1:1
    else
        d[1, 1] = 2
        d[2:n, 1] = 1 .+ view(dp, :, 1)
    end
    for i in 2:m
        if diags[i]
            d[1, i] = i
            d[2:n, i] .= 1 .+ view(dp, :, i - 1)
        else
            d[1, i] = 1 + d[1, i - 1]
            d[2:n, i] .= 1 .+ min.(view(dp, : i), view(d, 2:n, i - 1))
        end
    end
    d
end


    k = findfirst(diags)
    isnothing(k) && return [[(1 .+ dp) (2 .+ view(dp, :, m - 1))]; (2:m+1)']


function distances(diags::BitVector, d_prev)
    n = length(diags)
    d = [Matrix{Int}(undef, n + 1 - i, n) for i in 1:n]

    d[i][k + 1, j + 1] = 1 + min(dp[i][k, j + 1], d[i][k + 1, j])
    or
    d[i][k + 1, j + 1] = 1 + dp[i][k, j]

    d[i][1, j + 1] = 1 + d[i][j]
    or
    d[i][1, j + 1] = j + 1

    d[i][1, 1] = 2 - diags[i]

    d[i][k + 1, 1] = 1 + dp[i][k, 1]
    or
    d[i][k + 1, 1] = k + 1

    d[n][k + 1, 1] = d

    d[n] .= 2:n+1
    first_diag = findfirst(diags)
    if !isnothing(first_diag)
        d[n][first_diag:n] .-= 1
        d[n] .= 2:n+1
    else
        d[j][1:first_diag] .= 2:first_diag
    for k in 1:n
        if diags[n - k + 1]
            dn[k:n] .= k:n
            break
        end
        dn[k]
    d[n][

    for (i, di) in enumerate(d)
        d_above = 1:n
        for j in size(di)[1]
            dij = view(di, :, j)
            if diags[j]
                dij .= dprev
            else
                dij .= min(.d_above
            d[1, :] .= 1:n



function grow(s::Square)
    for d, v in zip(s.d, grow_indices(s)), (w, (a, b, c)) in placed_triples(s, v)
        if d[

    addable = Set(1:s.n+1)
    for (a, b, c) in triangles
        for ul in CartesianIndex(1, s.n - b + 2):CartesianIndex(s.n - a + 2, s.n - b + 2)
            remove blocked
            add disjunction
        end
    end
            pick smallest (?) set s.t. all disjunctions are satistifed

function bench(n)
    d = Vector{Matrix{Int}}()
    for i in 1:n
        d = distances(rand(Bool, i), d)
    end
    d
end

=#

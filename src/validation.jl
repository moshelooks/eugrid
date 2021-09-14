function sps(diags::AbstractMatrix{Bool}, edges::Bool=false)::Matrix{Int}
    d = Matrix{Int}(undef, size(diags) .+ 1)
    d[:, 1] = 0:size(diags)[1]
    d[1, :] = 0:size(diags)[2]
    for i in CartesianIndices(diags)
        d[i + onexy] = 1 + (diags[i] ? d[i] : min(d[i+onex], d[i+oney]))
    end
    edges && return d
    d[2:end, 2:end]
end

function iseugrid(diags::AbstractMatrix{Bool},
                  ts::Vector{Triple}=all_triples(maximum(size(diags))))::Bool
    for i in CartesianIndices(diags)
        d = sps(diags[i:CartesianIndex(size(diags))])
        x, y = size(d)
        for t in ts
            if t.a <= x && t.b <= y
                d[t.a, t.b] != t.c && return false
            end
        end
    end
    true
end

function maxbad(diags::AbstractMatrix{Bool},
                ts::Vector{Triple}=all_triples(maximum(size(diags))))
    mb = 0
    for i in CartesianIndices(diags)
        d = sps(diags[i:CartesianIndex(size(diags))])
        x, y = size(d)
        for t in ts
            if t.a <= x && t.b <= y
                mb = max(mb, abs(d[t.a, t.b] - t.c))
            end
        end
    end
    println(mb)
    mb
end


function fromint(n, x)::BitMatrix
    diags = BitMatrix(undef, n, n)
    for m in 0:n^2-1
        diags[m+1] = (x >> m) % 2 == 1
    end
    diags
end


brute(n) = (fromint(n, i-1) for i in Int128(1):Int128(2)^n^2)

function trace(kernel, basis, dmat, dbg=false)
    dset = Set(findall(dmat))
    o = Onion(kernel, basis, only(Set(size(dmat))))
    for ribbon in o.ribbons
        diags = intersect(dset, ribbon)
        dbg && println(diags)
        while true
            diags == popfirst!(o.solvers[end]) && break
        end
        didwrap = wrap!(o, diags)
        @assert didwrap
        if dbg && !iscomplete(o)
            a = only(o.solvers[end].stack)
            println(a.clauses)
            println(a.free)
            println(a.affirmed)
            println()
        end
    end
end

function trace_blocked(kernel, basis, dmat, dbg=false)
    dset = Set(findall(dmat))
    o = Onion(kernel, basis, only(Set(size(dmat))))
    for ribbon in o.ribbons
        diags = intersect(dset, ribbon)
        dbg && println(diags)
        while true
            #isempty(o.solvers[end]) && return o
            seen = popfirst!(o.solvers[end])
            diags == seen && break
            if length(o.solvers) < length(o.ribbons)
                ds = exterior_distances(o.membranes[end], seen)
                cs = constraints(o.box, ds)
                #println(cs)
                if !issatisfiable(cs)
                    #println("XXX ", violations(cs))
                    #=for v in violations(cs)
                        if v.u == Atom(4, 3)
                            return cs
                        end
                    end=#
                    bs = blockers(cs, o.membranes[end].interior, seen)
                    block!(o.solvers[end], bs)
                    for (b, v) in zip(bs, violations(cs))
                        if isblocked(Assignment([], Set(), diags), b)
                            println(b)
                            println(v)
                            println(seen)
                            return cs
                        end
                    end
                end
            end
        end
        didwrap = wrap!(o, diags)
        @assert didwrap
        if dbg && !iscomplete(o)
            a = only(o.solvers[end].stack)
            #println(a.clauses)
            #println(a.free)
            #println(a.affirmed)
            println()
        end
    end
end

function circ(diags, r=size(diags)[1]-2)
    x = sps(diags) .<= r
    [rot180(x) rotl90(x); rotr90(x) x]
end

function ecirc(n, r=n)
    x = BitMatrix(undef, n, n)
    for i in CartesianIndices(x)
        x[i] = sqrt(i[1]^2+i[2]^2) < r
    end
    x
end


function arc(diags, step, r=step)
    x = BitMatrix(undef, size(diags))
    for i in onexy:Atom(step, step):Atom(size(diags))
        br = min(Atom(size(diags)), i+Atom(step-1, step-1))
        x[i:br] .= sps(diags[i:br]) .< r
    end
    x
end

function earc(n, step)
    repeat(ecirc(step), outer=(div(n, step), div(n, step)))
end

function errstats(diags, step)
    es = Float64[]
    for i in onexy+Atom(3*step,3*step):Atom(size(diags) .- step)
        s = sps(diags[i:i+Atom(step, step)])
        for j in 1:step-1
            k = isqrt(step^2-j^2)
            push!(es, abs(s[j, k] - sqrt(j^2+k^2)))
        end
    end
    es
end

function geodesics(diags)::BitMatrix
    d = sps(diags, true)
    g = falses(size(d))
    stack = [Atom(size(d))]
    while !isempty(stack)
        i = pop!(stack)
        g[i] && continue
        g[i] = true
        minimum(i.I) == 1 && continue
        tl = i - onexy
        l = i - oney
        t = i - onex
        if diags[tl]
            push!(stack, tl)
            d[l] == d[tl] && push!(stack, l)
            d[t] == d[tl] && push!(stack, t)
        else
            d[l] <= d[t] && push!(stack, l)
            d[t] <= d[l] && push!(stack, t)
        end
    end
    g[1:findlast(g[:, 1]), 1] .= true
    g[1, 1:findlast(g[1,:])] .= true
    g
end

#=
struct Geodesic
    lambda::Int
    nodes::Set{Atom}
end


extend(g::Geodesic, a::Atom) = Geodesic(g.lambda + 1, union(g.nodes, [a]))

function nodiag(l::Geodesic, t::Geodesic, a::Atom)
    l.lambda < t.lambda && return extend(l, a)
    t.lambda < l.lambda && return extend(t, a)
    Geodesic(l.lambda + 1, union(l.nodes, t.nodes, [a]))
end

function geodesics(diags)::Matrix{Geodesic}
    ags = Matrix{Geodesic}(undef, size(diags) .+ 1)
    for i in 1:size(diags)[1]+1
        ags[i, 1] = Geodesic(i-1, Set(onexy:Atom(i, 1)))
    end
    for i in 2:size(diags)[2]+1
        ags[1, i] = Geodesic(i-1, Set(onexy:Atom(1, i)))
    end
    for i in onexy:Atom(size(diags))
        br = i + onexy
        ags[br] = diags[i] ? extend(ags[i], br) : nodiag(ags[i+onex], ags[i+oney], br)
    end
    ags
end

function ngeo(
    diags,
    ngeos::Vector{Vector{Int}}=[Vector{Int}() for _ in 1:sum(size(diags))])::
        Vector{Vector{Int}}
    for g in geodesics(diags)
        g.lambda > 0 && push!(ngeos[g.lambda], length(g.nodes))
    end
    ngeos
end

import Statistics

function all_ngeos(diags, k=maximum(sps(diags)))
    kmax = CartesianIndex(k, k)
    ngeos = [Vector{Int}() for _ in 1:sum(kmax.I)]
    for i in CartesianIndices(diags)
        ngeo(view(diags, i:min(i + kmax - onexy, CartesianIndex(size(diags)))), ngeos)
    end
    ngeos[1:k]
end
=#

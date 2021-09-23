function sps(diags::AbstractMatrix{Bool}, edges::Bool=false)::Matrix{Int}
    d = Matrix{Int}(undef, size(diags) .+ 1)
    d[:, 1] = 0:size(diags)[1]
    d[1, :] = 0:size(diags)[2]
    for i in CartesianIndices(diags)
        @inbounds d[i + onexy] = 1 + (diags[i] ? d[i] : min(d[i+onex], d[i+oney]))
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

#=
struct Geodesics
    x::CartesianIndex
    y::CartesianIndex
    d::Matrix{Int}
    g::BitMatrix
end
=#

function geodesics(diags)::BitMatrix
    d = sps(diags, true)
    g = falses(size(d))
    g[1, 1] = true
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

antidiag(diags, x::Atom) = Atom(2 + size(diags)[1] - x[1], x[2])

function distance(diags, x::Atom, y::Atom)::Int
    if y[2] < x[2]
        x, y = y, x
    end
    if y[1] < x[1]
        diags = view(diags, antidiag(diags, x):antidiag(diags, y)-onexy)
    else
        diags = view(diags, x:y-onexy)
    end
    sps(diags, true)[end]
end

function geodesics(diags, x::Atom, y::Atom)
    if y[2] < x[2]
        x, y = y, x
    end
    if y[1] < x[1]
        [Atom(x[1] - i[1] + 1, x[2] + i[2] - 1)
         for i in findall(geodesics(
             view(diags, antidiag(diags, x):antidiag(diags, y)-onexy)))]
    else
        Ref(x - onexy) .+ findall(geodesics(view(diags, x:y-onexy)))
    end
end

function circle(diags, center::Atom, radius::Int)::Set{Atom}
    pts = Set{Atom}([center+i*radius for i in (onex, oney, -onex, -oney)])
    br = Atom(size(diags))
    atrad(a, b, c) =
        findall(==(radius),
                sps(view(diags, min(max(a, onexy), br):b:min(max(c, onexy), br))))

    union!(pts, Ref(center) .+ atrad(center, onexy, center+Atom(radius-1, radius-1)))
    union!(pts, Ref(center) .- atrad(center-onexy, -onexy, center-Atom(radius,radius)))

    tcenter = antidiag(diags, center)
    union!(pts, [center + Atom(-x[1], x[2])
                 for x in atrad(tcenter, onexy, tcenter+Atom(radius-1, radius-1))])
    union!(pts, [center + Atom(x[1], -x[2])
                 for x in atrad(tcenter-onexy, -onexy, tcenter-Atom(radius,radius))])
    pts
end

function eqtriangle(diags, x::Atom, side::Int)
    pts = circle(diags, x, side)
    y = rand(pts)
    z = rand(intersect!(pts, circle(diags, y, side)))
    pts = Set{Atom}(geodesics(diags, x, y))
    union!(pts, geodesics(diags, y, z), geodesics(diags, z, x))
end

function samp(diags)
    x, y = rand(1:6145, 2)
    count(==(2048), sps(view(diags, x:x+2047, y:y+2047), true))
end

function sqrt3(diags, side::Int, x::Atom)
    pts = circle(diags, x, side)
    y = rand(pts)
    z = rand(intersect!(pts, circle(diags, y, side)))

    xm = Int(side/2)
    m = rand(intersect!(circle(diags, x, xm), circle(diags, y, xm)))
    distance(diags, m, z) / xm
end

sqrt3(diags, side::Int) =
    sqrt3(diags, side, rand(Atom(2*side+1, 2*side+1):Atom(size(diags) .- (2*side - 1))))

#=norm(x, y) = iisqrt(x^2 + y^2)
erf(x) = abs(x)

function choose_diag!(a::Atom, tl, l, t, br, avoid=nothing)::Bool
    n = length(l)
    br[1] = l[1] + 1
    br[n+1] = t[n] + 1
    br = view(br, 2:n)
    br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1

    !isnothing(avoid) && avoid[a] && return false

    s = 0.0

    a2sq = a[2]^2
    L = 0.0
    al = 0.0
    for i in 1:a[1]
        dsq = i^2 + a2sq
        de = sqrt(dsq)
        H = (abs(de - br[i]) - abs(de - tl[i] - 1)) / de
        ah = atan(i / a[2])

        delta = (H - L) / (ah - al)
        c = (L * ah - H * al) / (H - L)
        s += delta * (sqrt(c^2+2*c*ah+ah^2+1) * (c + ah) + asinh(c + ah) -
            (sqrt(c^2+2*c*al+al^2+1) * (c + al) + asinh(c + al)))
        #s += L * ah - H * al + 0.5 * (H - L) * (ah * sqrt(ah^2+1) + asinh(ah) - al * sqrt(al^2+1) - asinh(al))  / (ah - al)
        L = H
        al = ah
    end

    a1sq = a[1]^2
    for i in a[1]+1:a[1]+a[2]-1
        j = a[1]+a[2]-i
        dsq = a1sq + j^2
        de = sqrt(dsq)
        H = (abs(de - br[i]) - abs(de - tl[i] - 1)) / de
        ah = atan(a[1] / j)

        delta = (H - L) / (ah - al)
        c = (L * ah - H * al) / (H - L)
        s += delta * (sqrt(c^2+2*c*ah+ah^2+1) * (c + ah) + asinh(c + ah) -
            (sqrt(c^2+2*c*al+al^2+1) * (c + al) + asinh(c + al)))
        #s += L * ah - H * al + 0.5 * (H - L) * (ah * sqrt(ah^2+1) + asinh(ah) - al * sqrt(al^2+1) - asinh(al))  / (ah - al)
        #s += L * ah - H * al + 0.5 * (H - L) * (ah^2 - al^2) / (ah - al)
        L = H
        al = ah
    end

    H = 0
    ah = pi/2
    #s += L * ah - H * al + 0.5 * (H - L) * (ah^2 - al^2) / (ah - al)
    #s += L * ah - H * al + 0.5 * (H - L) * (ah * sqrt(ah^2+1) + asinh(ah) - al * sqrt(al^2+1) - asinh(al))  / (ah - al)
    delta = (H - L) / (ah - al)
    c = (L * ah - H * al) / (H - L)
    s += delta * (sqrt(c^2+2*c*ah+ah^2+1) * (c + ah) + asinh(c + ah) -
        (sqrt(c^2+2*c*al+al^2+1) * (c + al) + asinh(c + al)))

    if s >= 0
        br .= view(tl, 1:n-1) .+ 1
        true
    else
        false
    end
end=#

erf(x) = abs(x)

function norm(x, y, alpha=-2)
    #zx = exp(alpha * x)
    #zy = exp(alpha * y)
    #(x * zx + y*zy) / (zx + zy)
    #(Float64(x)^alpha + Float64(y)^alpha) ^ (1/alpha)
    #1/(1/x + 1/y)
    #min(x,y)/max(x,y)
    #(min(x,y)/max(x,y)) * 1e-2
    #min(x,y) / (x + y)
    #1/max(x,y)
    #x+y
    #min(x,y)*de
    #max(abs(x + y - de), abs(max(x,y) - de))
    #x * y / de
    #x += 1
    min(x, y)
end

function choose_diag!(a::Atom, tl, l, t, br, avoid=nothing)::Bool
    n = length(l)
    br[1] = l[1] + 1
    br[n+1] = t[n] + 1
    br = view(br, 2:n)
    br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1

    !isnothing(avoid) && avoid[a] && return false

    s = 0.0
    z = 0.0

    epsilon = 0 #(a[1] * a[2]) * 1e-10

    a2sq = a[2]^2
    for i in 1:a[1]-1
        dsq = i^2 + a2sq
        de = sqrt(dsq)
        w = atan(a[2] / (dsq - 0.25)) / norm(i, a[2])
        z += w
        ds = (erf(de - br[i]) - erf(de - tl[i] - 1))
        if ds == 0
            dwin = tl[i] + 1 < br[i]
            ndwin = l[i+1] != t[i]
            dwin && (ds += epsilon)
            ndwin && (ds -= epsilon)
        end
        s += ds * w
    end

    a1sq = a[1]^2

    dsq = a1sq + a2sq
    de = sqrt(dsq)
    w = atan((a[1]+a[2]-0.5) / (2 * dsq - a[1] - a[2])) / norm(a[1], a[2])
    z += w
    ds = (erf(de - br[a[1]]) - erf(de - tl[a[1]] - 1))
    if ds == 0
        dwin = tl[a[1]] + 1 < br[a[1]]
        ndwin = l[a[1]+1] != t[a[1]]
        dwin && (ds += epsilon)
        ndwin && (ds -= epsilon)
    end
    s += ds * w

    for i in a[1]+1:a[1]+a[2]-1
        j = a[1]+a[2]-i
        dsq = a1sq + j^2
        de = sqrt(dsq)
        w = atan(a[1] / (dsq - 0.25)) / norm(a[1], j)
        z += w
        ds = (erf(de - br[i]) - erf(de - tl[i] - 1))
        if ds == 0
            dwin = tl[i] + 1 < br[i]
            ndwin = l[i+1] != t[i]
            dwin && (ds += epsilon)
            ndwin && (ds -= epsilon)
        end
        s += ds * w
    end

    #push!(zs, z)

    if s >= 0
        br .= view(tl, 1:n-1) .+ 1
        true
    else
        false
    end
end

#=
function teqs(ds)
    s = 0
    for i in 1:length(ds)
        for j in 1:i-1
            if ds[j] + i - j > ds[i]
                s += 1
            end
        end
    end
    s
end

function choose_diag!(a::Atom, tl, l, t, br, avoid=nothing)::Bool
    n = length(l)
    br[1] = l[1] + 1
    br[n+1] = t[n] + 1
    br = view(br, 2:n)
    br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1

    !isnothing(avoid) && avoid[a] && return false

    s = 0
    ds = vcat(a[2], br[1:a[1]])
    s += teqs(ds)
    ds = vcat(a[2], tl[1:a[1]] .+ 1)
    s -= teqs(ds)
    ds = vcat(a[1], br[a[1] + a[2]-1:-1:a[1]])
    s += teqs(ds)
    ds = vcat(a[1], tl[a[1] + a[2]-1:-1:a[1]] .+ 1)
    s -= teqs(ds)

    if s <= 0
        br .= view(tl, 1:n-1) .+ 1
        true
    else
        false
    end
end
=#

function choose_children!(diags, a::Atom, grandparents, parents, children, avoid=nothing)
    Threads.@threads for j in 1:size(children, 2)
        tl = view(grandparents, :, j)
        l = view(parents, :, j)
        t = view(parents, :, j+1)
        br = view(children, :, j)
        a_j = a + CartesianIndex(-1, 1) * (j - 1)
        diags[a_j] = choose_diag!(a_j, tl, l, t, br, avoid)
    end
end

function children!(diags, a::Atom, grandparents, parents, children)
    Threads.@threads for j in 1:size(children, 2)
        tl = view(grandparents, :, j)
        l = view(parents, :, j)
        t = view(parents, :, j+1)
        br = view(children, :, j)
        a_j = a + CartesianIndex(-1, 1) * (j - 1)

        n = length(l)
        br[1] = l[1] + 1
        br[n+1] = t[n] + 1
        br = view(br, 2:n)
        if diags[a_j]
            br .= view(tl, 1:n-1) .+ 1
        else
            br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1
        end
    end
end

diags_buffer(n::Int) = Array{Int, 3}(undef, 2*n+1, n+2, 3)

function grow_base!(buffer)
    n = size(buffer, 2) - 2
    diags = BitMatrix(undef, n, n)

    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        choose_children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end

    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)
    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end

    diags
end

function grow_lower!(buffer, diags, avoid=nothing)
    n = checksquare(diags)
    grandparents = view(buffer, 1:1, 1:1, 1) .= 0
    parents = view(buffer, 1:2, 1:2, 2) .= [0 1; 1 0]
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)

    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children, avoid)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end
end

function grow_diags(n::Int)::BitMatrix
    buffer = diags_buffer(n+1)
    base = grow_base!(buffer)
    diags = base[(n+1)*onexy:-onexy:2*onexy]
    grow_lower!(buffer, diags)
    diags
end

struct Grid
    diags::BitMatrix
    antidiags::BitMatrix
end

Base.getproperty(g::Grid, s::Symbol) = s === :n ? size(g.diags, 1) : getfield(g, s)

function grow_grid(n::Int)::Grid
    buffer = diags_buffer(2n)
    base = grow_base!(buffer)
    diags = base[(n+1)*onexy:2n*onexy]

    avoid = falses(2n, 2n)
    avoid[(n+1)*onexy:2n*onexy] .= view(diags, n:-1:1, :)
    grow_lower!(buffer, base, avoid)
    antidiags = base[(n+1)*onexy:2n*onexy]

    Grid(diags, antidiags)
end

function score(g)
    g = BitMatrix(g)
    n = size(g)[1]
    step = Int(n / 8)
    sum(arc(g, step) .!= earc(n, step)) / n^2
end

function grow_simple(n::Int, leq=false)::BitMatrix
    ds = Matrix{Int}(undef, n+1, n+1)
    ds[:, 1] = ds[1, :] = 0:n
    diags = BitMatrix(undef, n, n)
    for i in CartesianIndices(diags)
        tl = ds[i]
        l = ds[i+onex]
        t = ds[i+oney]
        de = sqrt(i[1]^2 + i[2]^2)
        d_diag = tl + 1
        d_nodiag = min(l, t) + 1
        if (abs(de - d_diag) < abs(de - d_nodiag) ||
            (leq && abs(de - d_diag) == abs(de - d_nodiag)))
            ds[i + onexy] = d_diag
            diags[i] = true
        else
            ds[i + onexy] = d_nodiag
            diags[i] = false
        end
    end
    diags
end

function minscore(g, n, skip=1)
    to = Atom(n-1,n-1)
    ix = onexy:Atom(skip, skip):Atom(size(g))-to
    scores = [score(g[i:i+to]) for i in ix]
    s, i = findmin(scores)
    s, ix[i]
end

function arcs(diags, center, radii, w=1)
    radii = Set{Int}(Iterators.flatten(r .+ (0:w-1) for r in radii))
    w = maximum(radii)
    d = sps(diags[center:center+(w-1)*onexy])
    .![d[i] in radii for i in Atoms(w)][w:-1:1, :]
end

using Statistics

randgrid(n, p) = Grid(rand(Float64, (n, n)) .< p, rand(Float64, (n, n)) .< p)

function gruntle(n, k)
    ps = 0.18:0.001:0.2
    ps[argmin([median([abs(sqrt(2)*n-sps(rand(Float64, (n, n)) .< p)[end]) for _ in 1:k])
               for p in ps])]
end


function gruntle2(n, k)
    ps = 0.0:0.005:0.19
    mangles = [rmangle(randgrid(n, p), k)[2] for p in ps]
    x, ix = findmin(mangles)
    x, ps[ix]
end

function ngs(g, k)
    a = Analysis(g)
    for _ in 1:k
        sample!(a)
    end
    xs = [x.second/x.first for x in a.ngs]
    alpha = sum(log2.(last.(a.ngs))) / sum(-1 .+ log2.(first.(a.ngs)))
    alpha, minimum(xs), mean(xs), maximum(xs)
end

function mungs(g, k, a=rand(Atoms(513)))
    g = Grid(g.diags[a:a+onexy*511], g.antidiags[514-a[1]:1025-a[1],a[2]:a[2]+511])
    @assert maximum(g.diags .+ g.antidiags[512:-1:1,:]) == 1
    ngs(g, k), a
end

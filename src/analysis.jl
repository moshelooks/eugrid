struct Analysis
    g::Grid
    seed::Int
    rng::StableRNG
    es::Vector{Int}
    ngs::Vector{Pair{Int, Int}}
    s3s::Vector{Pair{Int, Float64}}

    Analysis(g::Grid, seed::Int=1) =
        new(g, seed, StableRNG(seed), Int[], Pair{Int, Int}[], Pair{Int, Float64}[])
end


function tests(g::Grid; seed=1, k=100)
    rng = StableRNG(seed)
    vs = rand(vertices(g), k)
end

function sample!(a::Analysis)
    atoms = Atoms(a.g.n)
    x = rand(a.rng, atoms)
    dx = sps(a.g, x)
    e = maximum(dx)
    push!(a.es, e)

    y = x
    while dx[y] < 4
        y = rand(a.rng, atoms)
    end
    side = dx[y]
    dy = sps(a.g, y, side)

    #println(x)
    #println(y)

    ngeo = 0
    for i in atoms
        dx[i] + dy[i] == side && (ngeo += 1)
    end

    #println(side=>ngeo)

    push!(a.ngs, side=>ngeo)
    #=
    if isodd(side) || 3*side > a.g.n
        while dx[y] < 4 || isodd(dx[y]) || 3*dx[y] > a.g.n
            y = rand(a.rng, atoms)
        end
        side = dx[y]
        dy = sps(a.g, y, side)
    end

    xm = Int(side / 2)
    m = rand(a.rng, [i for i in atoms if dx[i] == dy[i] == xm])

    if isempty([i for i in atoms if dx[i] == dy[i] == side])
        println(x)
        println(y)
    end

    z = rand(a.rng, [i for i in atoms if dx[i] == dy[i] == side])
    zm = sps(a.g, z, side)[m]
    s3 = zm / xm

    if s3 > 10
        println(x)
        println(y)
        println(z)
    end


    push!(a.s3s, side=>s3)
    =#
    a
end



#=

shift(x, a::Atom) = CircularArray(x)[a:a+Atom(size(x) .- 1)].data
unshift(x, a::Atom) = rotate(x, Atom(2, 2) - a)

antishift(x, a::Atom) = shift(x, Atom(2 - a[1], a[2]))
antiunshift(x, a::Atom) =
    CircularArray(x)[a[1]:-1:a[1]-size(x,1)+1,2-a[2]:1-a[2]+size(x,2)].data

dsps(diags, a::Atom) =
    unshift(sps(view(shift(diags, a), onexy:Atom(size(diags) .- 1)), true), a)

adsps(diags, a::Atom) =
    antiunshift(sps(view(antishift(diags, a), onexy:Atom(size(diags) .- 1)), true), a)

function sps(g::Grid, a::Atom, m=div(2*g.n, 3, RoundUp))::Matrix{Int}
    n = checksquare(g.diags)
    d = CircularArray(fill(2n, size(g.diags)))

    cd = CircularArray(g.diags)
    dv = view(d, a:a+onexy*m)
    dv .= min.(dv, sps(view(cd, a:a+onexy*(m-1)), true))

    dv = view(d, a:-onexy:a-onexy*m)
    dv .= min.(dv, sps(view(cd, a-onexy:-onexy:a-onexy*m), true))

    cd = CircularArray(g.antidiags)
    step = oney - onex
    dv = view(d, a+step:step:a+step*m)
    dv .= min.(dv, sps(view(cd, 2-a[1]:1-a[1]+m, a[2]:a[2]+m-1)))

    step = onex - oney
    dv = view(d, a+step:step:a+step*m)
    dv .= min.(dv, sps(view(cd, 1-a[1]:-1:2-a[1]-m, a[2]-1:-1:a[2]-m)))

    d.data
end

function n_geodesics(g::Grid, x::Atom, y::Atom, dx=sps(g, x), dy=sps(g, y, dx[y]))
    side = dx[y]
    ngeo = 0
    for i in Atoms(g.n)
        dx[i] + dy[i] == side && (ngeo += 1)
    end
    ngeo
end

function max_geodesics(g::Grid)
    spss = sps.(Ref(g), Atoms(g.n))
    ngss = [n_geodesics(g, x, y, spss[x], spss[y]) for x in Atoms(g.n), y in Atoms(g.n)]
    maximum(ngss)
end

function geo_alpha(g::Grid)
    spss = sps.(Ref(g), Atoms(g.n))
    num = denom = 0.0
    for x in Atoms(g.n), y in Atoms(g.n)
        side = spss[x][y]
        side < 4 && continue
        num += log2(n_geodesics(g, x, y, spss[x], spss[y]))
        denom += log2(side) - 1
    end
    num / denom
end

function sample!(a::Analysis)
    atoms = Atoms(a.g.n)
    x = rand(a.rng, atoms)
    dx = sps(a.g, x)
    e = maximum(dx)
    push!(a.es, e)

    y = x
    while dx[y] < 4
        y = rand(a.rng, atoms)
    end
    side = dx[y]
    dy = sps(a.g, y, side)

    #println(x)
    #println(y)

    ngeo = 0
    for i in atoms
        dx[i] + dy[i] == side && (ngeo += 1)
    end

    #println(side=>ngeo)

    push!(a.ngs, side=>ngeo)
    #=
    if isodd(side) || 3*side > a.g.n
        while dx[y] < 4 || isodd(dx[y]) || 3*dx[y] > a.g.n
            y = rand(a.rng, atoms)
        end
        side = dx[y]
        dy = sps(a.g, y, side)
    end

    xm = Int(side / 2)
    m = rand(a.rng, [i for i in atoms if dx[i] == dy[i] == xm])

    if isempty([i for i in atoms if dx[i] == dy[i] == side])
        println(x)
        println(y)
    end

    z = rand(a.rng, [i for i in atoms if dx[i] == dy[i] == side])
    zm = sps(a.g, z, side)[m]
    s3 = zm / xm

    if s3 > 10
        println(x)
        println(y)
        println(z)
    end


    push!(a.s3s, side=>s3)
    =#
    a
end


#=
function eccentricity(g::Grid, x::Atom)::Float64

end


function eccentricity(rng::StableRNG, g::Grid, N::Int)::Vector{Float64}

end





    for i in a:a+onexy*(m-1)
    d
end

                         1:end-1, 1:end-1), a)
    ad = sps(CircularArray(g.antidiags)[a[1]:-1:a[1]-g.n+1, a[2]:a[2]+g.n-1])[
    min.(d, ad)



function distance(g::Grid, x::Atom, y::Atom)::Int
    if y[2] < x[2]
        x, y = y, x
    end
    sps(y[1] < x[1] ?
        view(g.antidiags, antidiag(diags, x):antidiag(diags, y)-onexy) :
        view(g.diags, x:y-onexy),
        true)[end]
end

clamp(a::Atom, g::Grid) = min(max(a, onexy), Atom(size(a.diags)))
#=
function quarter(g::Grid, center::Atom, radius::Int, lower::Bool, right::Bool, open::Bool)
    if lower
        if right
            diags =
    dx = dy = radius
    upper && (dx *= )
    [center + Atom(i.I .* rot) for i in findall(==(radius), sps(diags, closed))]
end
=#


function circle(g::Grid, center::Atom, radius::Int)::Set{Atom}
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

function ngeodesics(rng::StableRNG, g::Grid, N::Int)::Vector{Pair{Int, Int}}

end

function sqrt3(rng::StableRNG, g::Grid, side::Int, N::Int)::Vector{Float64}

end

function analyze(g::Grid, N::Int seed::Int=1)
    rng = StableRNG(seed)
    es = eccentricity(StableRNG(abs(rand(rng, Int))), g, N)
    ngs = ngeodesics(StableRNG(abs(rand(rng, Int))), g, N)
    s3s = sqrt3(StableRNG(abs(rand(rng, Int))), g, Int(round(Statistics.mean(
    Analysis(seed
=#

function dmin(a, b, n)
    a, b = minmax(a, b)
    min(b - a, a + n - b)
end


function mangle_stats(g::Grid, a::Atom)
    m = Int(g.n / 2) - 1
    d = [Vector{Int}() for _ in 1:4]
    for i in findall(sps(g, a, m) .== m)
        i -= a
        if i[1] > m
            i = Atom(i[1] - g.n, i[2])
        elseif i[1] < -m
            i = Atom(i[1] + g.n, i[2])
        end
        if i[2] > m
            i = Atom(i[1], i[2] - g.n)
        elseif i[2] < -m
            i = Atom(i[1], i[2] + g.n)
        end
        if i[1] == m
            push!(d[1], i[2])
        elseif i[1] == -m
            push!(d[2], i[2])
        end
        if i[2] == m
            push!(d[3], i[1])
        elseif i[2] == -m
            push!(d[4], i[1])
        end
    end
    d
end

function mangle(g::Grid, a::Atom)
    m = Int(g.n / 2) - 1
    d = 0
    for (lo, hi) in extrema.(mangle_stats(g, a))
        d = max(d, min(atand((-lo-1)/m) + atand(hi/m), atand(-lo/m) + atand((hi+1)/m)))
    end
    d
end
#=


    d = 0


    d = 0
    for i in findall(sps(g, a, m) .== m)
        i -= a
        if i[1] > m
            i = Atom(i[1] - m, i[2])
        end
        if i[2] > m
            i = Atom(i[1], i[2] -1)
        end
        if i[1] == m
            push!(d[1], i[2])
        elseif i[1] == -m
            push!(d[2], i[2])
        elseif i[2] == m
            push!(d[3], i[1])
        elseif i[2] == -m
            push!(d[4], i[1])
        end
    end
    d


        if a[1] - i[1] == m
            if


        x, y = minmax(dmin.(a.I, i.I, Ref(g.n))...)
        @assert y <= m
        y != m && continue

            ix = 1
        elseif abs(a[2] - i[2]) == m
            ix = 2
        elseif dmin(a,



        d = max(d, x)
    end
    atand((d + 1) / m)
end
=#
using Statistics

function rmangle(g::Grid, k, seed=1)
    rng = StableRNG(seed)
    xs = [mangle(g, rand(rng, Atoms(g.n))) for _ in 1:k]
    minimum(xs), mean(xs), maximum(xs)
end

function circs(g::Grid, n)
    w = Int(g.n / n)
    r = div(w, 3)
    x = trues(g.n, g.n)
    for i in onexy:w*onexy:onexy*g.n
        x[findall(sps(g, i+Int(w/2)*onexy, r) .< r)] .= false
    end
    x
end

using DataFrames, GLM

function geocount(g::Grid, nsamples=100, seed=1)
    rng = StableRNG(seed)
    mingeo = 4
    maxgeo = 2 * Int(floor(0.5 * 0.5 * 0.7 * g.n))
    atoms = Atoms(g.n)

    X = Float64[]
    Y = Float64[]

    for m in mingeo:2:maxgeo
        ngeo_sum = ngeo_sum2 = ngeo_logsum = ngeo_logsum2 = 0.0

        for _ in 1:nsamples
            x = rand(rng, atoms)
            dx = sps(g, x, m)
            y = rand(rng, findall(==(m), dx))

            ngeo = n_geodesics(g, x, y, dx)
            ngeo_sum += ngeo
            ngeo_sum2 += ngeo^2

            wk = log2(ngeo)
            ngeo_logsum += wk
            ngeo_logsum2 += wk^2
        end

        av = ngeo_sum / nsamples
        av2 = ngeo_sum2 / nsamples
        stdev = sqrt(av2 - av*av)
        #println("$m $av $stdev")

        av = ngeo_logsum/ nsamples
        av2 = ngeo_logsum2 / nsamples
        stdev = sqrt(av2 - av*av)
        println("$(log2(m)) $av $stdev")
        push!(X, log2(m))
        push!(Y, av)
    end
    data = DataFrame(X=X, Y=Y)
    ols = lm(@formula(Y ~ X), data)
    coef(ols)[2]
end
=#

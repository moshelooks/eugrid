"""
The basic idea here is to separate out the generation of random variates from the
computational heavy lifting so that the latter can be multi-threaded and we'll still have
reproducible results regardlessof how many threads there are.

10,000 xs for eccentricity

~ 0.35 * n * 100 ys for geocount

"""

struct Sample
    seed::Int

    x::Vertex
    y::Vertex
    z::Vertex

    eccentricity_x::Int
    euclidean_eccentricity_x::Float64

    dxy::Int
    geocount_xy::Int

    c::Vertex
    dcz::Int
end

function Sample(g::Grid, seed::Int, dxy::Int)::Sample
    rng = StableRNG(seed)

    x = rand(rng, vertices(g)[onexy*(dxy+1):onexy*(g.n-dxy)])
    dx = sps(g, x)

    eccentricity_x = eccentricity(g, x, dx)
    euclidean_eccentricity_x = euclidean_eccentricity(g.n, x)

    dxy == 0 && return Sample(
        seed, x, x, x, eccentricity_x, euclidean_eccentricity_x, dxy, 1, x, 0)

    y = rand(rng, circle_points(g, x, dxy, dx))
    @assert y != x

    dy = sps(g, y, dxy + 1)
    geocount_xy = sum(geodesics(g, x, y, dx, dy))

    dxy < (g.n - 1) / 16 && return Sample(
        seed, x, y, y, eccentricity_x, euclidean_eccentricity_x, dxy, geocount_xy, y, 0)

    z = rand(rng, two_circle_points(g, x, y, dxy, dx, dy, strict=false))
    @assert z != y

    c = rand(rng, midpoints(g, x, y, dx, dy))
    dcz = d2i(distance(g, c, z), RoundDown)

    Sample(seed, x, y, z, eccentricity_x, euclidean_eccentricity_x, dxy, geocount_xy, c, dcz)
end

delta_eccentricity(samples)::Vector{Float64} =
    [s.eccentricity_x - s.euclidean_eccentricity_x for s in samples]

function geodesic_model(samples)
    samples = [s for s in samples if s.x != s.y]
    xs = [log2(s.dxy) for s in samples]
    ys = [log2(s.geocount_xy) for s in samples]
    GLM.lm(@GLM.formula(Y ~ X), DataFrames.DataFrame(X=xs, Y=ys))
end

function geodesic_exponent(samples)::Float64
    GLM.coef(geodesic_model(samples))[2]
end

delta_sqrt3(samples)::Vector{Float64} =
    [s.dcz / s.dxy - sqrt(3) / 2 for s in samples if s.y != s.z]

function score(samples)
    de = delta_eccentricity(samples)
    mu_de = Statistics.mean(de)
    ds = delta_sqrt3(samples)
    [mu_de, mu_de + Statistics.std(de),
     geodesic_exponent(samples),
     isempty(ds) ? NaN : Statistics.mean(ds)]
end

function sample(g::Grid, k::Int; seed::Int=1, min_geo::Int=4, max_geo::Int=2*div(7*g.n, 40))
    @assert max_geo >= min_geo
    nzeros = max(0, round(Int, k * (100 / (max_geo - min_geo + 1) - 1), RoundUp))
    dxys = repeat([zeros(Int, nzeros); min_geo:max_geo], k)
    samples = Vector{Sample}(undef, length(dxys))
    Threads.@threads for i in eachindex(dxys, samples)
        samples[i] = Sample(g, i + (seed - 1) * length(dxys), dxys[i])
    end
    samples
end

function H(g::Grid, k)
    denom = 1.0 #expected_euclidean_eccentricity(g.n)
    samples = [(s.eccentricity_x - s.euclidean_eccentricity_x) / denom
               for s in sample(g, k, min_geo=0, max_geo=0)]
    Statistics.mean(samples), Statistics.std(samples)
end

function hmin(n::Int, rng=StableRNG(1))
    seed = rand(rng, 1:2^31)
    res = Optim.optimize(p->(H(randgrid(StableRNG(seed), n, p), 1)[1]^2), 0.0, 1.0)
    println(res)
    randgrid(StableRNG(seed), n, res.minimizer)
end


function randh(n::Int, k, rng=StableRNG(1))
    seed = rand(rng, 1:2^31)
    res = Optim.optimize(p->(H(randgrid(StableRNG(seed), n+1, p), 1)[1]^2), 0.0, 1.0)
    println(res)
    H(randgrid(StableRNG(seed), n+1, res.minimizer), k)
end

function sample_HG_AG(g::Grid, rng_base::StableRNG, k::Int)
    seeds = rand(rng_base, 1:2^32, k)
    hg = Vector{Float64}(undef, k)
    ag = Vector{Float64}(undef, k)
    vs = vertices(g)
    Threads.@threads for i in 1:k
        rng = StableRNG(seeds[i])
        v = rand(rng, vs)
        dv = sps(g, v)
        hg[i] = HG(g, v, dv)
        ag[i] = AG(g, v, dv)
    end
    hg, ag
end

function sample_NG_TG(g::Grid, d::Int, rng_base::StableRNG, k::Int; skip_tg=false)
    seeds = rand(rng_base, 1:2^32, k)
    ng = Vector{Int64}(undef, k)
    tg = Vertex{Float64}(undef, skip_tg ? 0 : k)
    vs = vertices(g)[onexy*(d+1):onexy*(g.n-d)]
    Threads.@threads for i in 1:k
        rng = StableRNG(seeds[i])

        x = rand(rng, vs)
        dx = sps(g, x)

        y = rand(rng, circle_points(g, x, dxy, dx))
        @assert y != x

        dy = sps(g, y, dxy + 1)
        ng[i] = NG(g, x, y, dx, dy)

        skip_tg && continue

        z = rand(rng, two_circle_points(g, x, y, dxy, dx, dy, strict=false))
        @assert z != y

        c = rand(rng, midpoints(g, x, y, dx, dy))
        dcz = d2i(distance(g, c, z), RoundDown)

        tg[i] = TG(g, d, dcz)
    end
    ng, tg
end

struct Sampling
    hg::Dict{Int, Vector{Float64}}
    ag::Dict{Int, Vector{Float64}}
    ng::Dict{Int, Dict{Int, Vector{Int}}}
    tg::Dict{Int, Dict{Int, Vector{Float64}}}
end

Sampling() = Sampling(Dict(), Dict(), Dict(), Dict())

function sample!(s::Sampling, g::Grid, rng::StableRNG, k::Int)
    hg, ag = sample_HG_AG(g, rng, min(k, 100)^2)
    append!(get!(()->Float64[], s.hg, g.n), hg)
    append!(get!(()->Float64[], s.ag, g.n), ag)

    ngs = get!(()->Dict{Int, Vector{Int}}(), s.ng, g.n)
    tgs = get!(()->Dict{Int, Vector{Float64}}(), s.tg, g.n)

    mingeo = 4
    mingeo_TG = 2*div(7*g.n, 200)
    @assert mingeo_TG >= mingeo
    maxgeo = 2*div(7*g.n, 40)

    for d in mingeo:2:maxgeo
        ng, tg = sample_NG_TG(g, d, rnd, k, d < mingeo_TG)
        append!(get!(()->Int[], ngs, d), ng)
        append!(get!(()->Float64[], tgs, d), tg)
    end
end

function rand_sampling(seed=1, nreplicates=25, k=100, ns=[2^i+1 for i in 6:10])
    rng = StableRNG(seed)
    s = Sampling()
    for n in ns
        for _ in 1:nreplicates
            g = randmin(n, rng)
            sample!(s, g, rng, k)
        end
    end
    s
end
#=
function diags_sampling(diags, seed=2, k=500, ns=[2^i+1 for i in 6:10])
    rng = StableRNG(seed)
    s = Sampling()
    for n in ns


            hg, ag = sample_HG_AG(g, rng, k^2)
            append!(ss.hg[n], hg)
            append!(ss.ag[n], ag)
            ng, tg

=#

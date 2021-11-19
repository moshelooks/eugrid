function sample_HG(rng_base::StableRNG, g::Grid, nsamples::Int)
    seeds = rand(rng_base, 1:2^32, nsamples)
    hg = Vector{Float64}(undef, nsamples)
    vs = vertices(g)
    Threads.@threads for i in 1:nsamples
        rng = StableRNG(seeds[i])
        v = rand(rng, vs)
        dv = sps(g, v)
        hg[i] = HG(g, v, dv)
    end
    hg
end

function sample_AG_NG_TG(rng_base::StableRNG, g::Grid, dxy::Int, nsamples::Int)
    seeds = rand(rng_base, 1:2^32, nsamples)
    ng_only = dxy < 2*div(7*g.n, 200)
    ag = Vector{Float64}(undef, ng_only ? 0 : nsamples)
    ng = Vector{Int64}(undef, nsamples)
    tg = Vector{Float64}(undef, ng_only ? 0 : nsamples)
    Threads.@threads for i in 1:nsamples
        rng = StableRNG(seeds[i])

        x = Vertex(rand(rng, 1+dxy:g.n-dxy), rand(rng, 1+dxy:g.n-dxy))
        dx = sps(g, x, dxy + 1)

        y = rand(rng, circle_points(g, x, dxy, dx))
        @assert y != x

        dy = sps(g, y, dxy + 1)
        ng[i] = NG(g, x, y, dx, dy)

        ng_only && continue

        ag[i] = AG(rng, g, dxy, x, dx)

        z = rand(rng, two_circle_points(g, x, y, dxy, dx, dy, strict=false))
        @assert z != y

        c = rand(rng, midpoints(g, x, y, dx, dy))
        dcz = d2i(distance(g, c, z), RoundDown)

        tg[i] = TG(g, dxy, dcz)
    end
    ag, ng, tg
end

struct Sampling
    hg::Dict{Int, Vector{Float64}}
    ag::Dict{Int, Dict{Int, Vector{Float64}}}
    ng::Dict{Int, Dict{Int, Vector{Int}}}
    tg::Dict{Int, Dict{Int, Vector{Float64}}}
end

Sampling() = Sampling(Dict(), Dict(), Dict(), Dict())

function geodesic_model(ngs::Dict{Int, Vector{Int}})
    xs = Float64[]
    ys = Float64[]
    for (d, ng_d) in ngs
        log2d = log2(d)
        for ng in ng_d
            push!(xs, log2d)
            push!(ys, log2(ng))
        end
    end
    GLM.lm(@GLM.formula(Y ~ X), DataFrames.DataFrame(X=xs, Y=ys))
end

function sample!(rng::StableRNG, s::Sampling, g::Grid, nsamples1::Int, nsamples2::Int)
    println(g.n, " hg")
    hg = sample_HG(rng, g, nsamples1)
    append!(get!(()->Float64[], s.hg, g.n), hg)

    ags = get!(()->Dict{Int, Vector{Float64}}(), s.ag, g.n)
    ngs = get!(()->Dict{Int, Vector{Int}}(), s.ng, g.n)
    tgs = get!(()->Dict{Int, Vector{Float64}}(), s.tg, g.n)

    mingeo = 4
    maxgeo = 2*div(7*g.n, 40)

    for d in mingeo:2:maxgeo
        println(g.n, " ag_ng_tg: ", d, " / ", maxgeo)

        ag, ng, tg = sample_AG_NG_TG(rng, g, d, nsamples2)
        append!(get!(()->Float64[], ags, d), ag)
        append!(get!(()->Int[], ngs, d), ng)
        append!(get!(()->Float64[], tgs, d), tg)
    end
end

function diags_sampling(diags, nsamples1=10000, nsamples2=1000, ns=[2^i+1 for i in 6:13];
                        seed=1)
    rng = StableRNG(seed)
    s = Sampling()
    for n in ns
        g = Grid(diags, n)
        sample!(rng, s, g, nsamples1, nsamples2)
    end
    s
end

randgrid(rng::StableRNG, n::Int, p::Float64) =
    Grid(rand(rng, Float64, (n-1, n-1)) .< p, rand(rng, Float64, (n-1, n-1)) .< p)

function randmin(rng_base::StableRNG, n::Int)
    seed = rand(rng_base, 1:2^32)
    Optim.optimize(0.0, 1.0) do p
        rng = StableRNG(seed)
        g = randgrid(rng, n, p)
        hg = sample_HG(rng, g, 128)
        Statistics.mean(hg .^ 2)
    end
end

function rand_sampling(nreplicates=25, nsamples1=10000, nsamples2=100,
                       ns=[2^i+1 for i in 6:11]; seed=2)
    rng = StableRNG(seed)
    s = Sampling()
    for n in ns
        res = randmin(rng, n)
        println(res)
        p = res.minimizer
        for _ in 1:nreplicates
            g = randgrid(rng, n, p)
            sample!(rng, s, g, nsamples1, nsamples2)
        end
    end
    s
end

function draw_hg(samples)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, subsample) in samples
        N = n^2
        x = log2(N)
        hs = delta_eccentricity(subsample) ./ x
        y = Statistics.mean(hs)
        yerr = Statistics.std(hs)
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!("\$\$\\log_2{N}\$\$")
    Plots.ylabel!("\$\$H_G\$\$")
end

function draw_ngeo(samples)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, subsample) in samples
        N = n^2
        x = log2(N)
        model = geodesic_model(subsample)
        y = GLM.coef(model)[2]
        yerr = GLM.stderror(model)[2]
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!("\$\$\\log_2{N}\$\$")
    Plots.ylabel!("\$\$\\gamma\$\$")
end

function draw_tg(samples)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, subsample) in samples
        N = n^2
        x = log2(N)
        ts = delta_sqrt3(subsample) ./ x
        y = Statistics.mean(ts)
        yerr = Statistics.std(ts)
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!("\$\$\\log_2{N}\$\$")
    Plots.ylabel!("\$\$T_G\$\$")
end

function randmin(n::Int, rng::StableRNG)
    seed = rand(rng, 1:2^32)
    res = Optim.optimize(0.0, 1.0) do p
        subrng = StableRNG(seed)
        g = randgrid(subrng, n, p)
        samples = sample(g, 1, seed=rand(subrng, 1:2^32), min_geo=0, max_geo=0)
        Statistics.mean(delta_eccentricity(samples))^2
    end
    println(res)
    subrng = StableRNG(seed)
    p = res.minimizer
    randgrid(subrng, n, p)
end

function sample_rand(seed=1, nreplicates=5, ns=[2^i+1 for i in 6:10], k=20)
    rng = StableRNG(seed)
    samples = Dict{Int, Vector{Sample}}(n=>Sample[] for n in ns)
    for n in ns
        for _ in 1:nreplicates
            g = randmin(n, rng)
            append!(samples[n], sample(g, k, seed=rand(rng, 1:2^32)))
        end
    end
    samples
end

function draw_rand(samples)
    draw_hg(samples)
    Plots.savefig("paper/images/rand_hg.svg")
    draw_ngeo(samples)
    Plots.savefig("paper/images/rand_ngeo.svg")
    draw_tg(samples)
    Plots.savefig("paper/images/rand_tg.svg")
end

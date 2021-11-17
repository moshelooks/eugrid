function draw_hg(s::Sampling)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, hg) in s.hg
        N = n^2
        x = log2(N)
        y = Statistics.mean(hg)
        yerr = Statistics.std(hg)
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!("\$\$\\log_2{N}\$\$")
    Plots.ylabel!("\$\$H_G\$\$")
end

function draw_gg(s::Sampling)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, ngs) in s.ng
        N = n^2
        x = log2(N)
        model = geodesic_model(ngs)
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

function draw_tg(s::Sampling)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, tgs) in s.tg
        N = n^2
        x = log2(N)
        ts = collect(Iterators.flatten(values(tgs)))
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

function draw_ag(s::Sampling)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, ag) in s.ag
        N = n^2
        x = log2(N)
        y = Statistics.mean(ag)
        yerr = Statistics.std(ag)
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!("\$\$\\log_2{N}\$\$")
    Plots.ylabel!("\$\$A_G\$\$")
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
    return res.minimizer
    subrng = StableRNG(seed)
    p = res.minimizer
    randgrid(subrng, n, p)
end

function draw_rand(s::Sampling)
    draw_hg(s)
    Plots.savefig("paper/images/rand_hg.svg")
    draw_gg(s)
    Plots.savefig("paper/images/rand_gg.svg")
    draw_tg(s)
    Plots.savefig("paper/images/rand_tg.svg")
    draw_ag(s)
    Plots.savefig("paper/images/rand_ag.svg")
end

function draw_gamma(s::Sampling)
    draw_hg(s)
    Plots.savefig("paper/images/gamma_hg.svg")
    #draw_gg(s)
    #Plots.savefig("paper/images/gamma_gg.svg")
    #draw_tg(s)
    #Plots.savefig("paper/images/gamma_tg.svg")
    draw_ag(s)
    Plots.savefig("paper/images/gamma_ag.svg")
end

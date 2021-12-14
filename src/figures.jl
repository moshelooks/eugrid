function draw_hg(s::Sampling)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, hg) in s.hg
        N = n^2
        x = log2(N)
        hg = (hg * x) / sqrt(n)
        y = Statistics.mean(hg)
        yerr = Statistics.std(hg)
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!(L"\log_2{N}")
    Plots.ylabel!(L"H_G")
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
    Plots.xlabel!(L"\log_2{N}")
    Plots.ylabel!(L"\gamma")
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
    Plots.xlabel!(L"\log_2{N}")
    Plots.ylabel!(L"T_G")
end

function draw_ag(s::Sampling)
    xs = Float64[]
    ys = Float64[]
    yerrs = Float64[]
    for (n, ags) in s.ag
        N = n^2
        x = log2(N)
        ag = collect(Iterators.flatten(values(ags)))
        y = Statistics.mean(ag)
        yerr = Statistics.std(ag)
        push!(xs, x)
        push!(ys, y)
        push!(yerrs, yerr)
    end
    Plots.scatter(xs, ys, yerror=yerrs, grid=false, legend=false, mc=:black, shape=:hline,
                  bg=:transparent, fg=:black, yguidefontsize=14, xguidefontsize=14)
    Plots.xlabel!(L"\log_2{N}")
    Plots.ylabel!(L"A_G")
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
    draw_gg(s)
    Plots.savefig("paper/images/gamma_gg.svg")
    draw_tg(s)
    Plots.savefig("paper/images/gamma_tg.svg")
    draw_ag(s)
    Plots.savefig("paper/images/gamma_ag.svg")
end

function draw_gamma_disordered(s::Sampling)
    draw_gg(s)
    Plots.savefig("paper/images/gamma_disordered_gg.svg")
end

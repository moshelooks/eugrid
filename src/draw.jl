name(i::CartesianIndex{2}) = "($(i[1])x$(i[2]))"

node(i::CartesianIndex{2}) =
    "  \\node $(name(i)) at ($(i[2]-1), $(i[1]-1)) [inner sep=0pt,minimum size=0pt,draw] {};"

nostyle(i, j) = ""

edge(i::CartesianIndex{2}, j::CartesianIndex{2}, style=nostyle) =
    "  \\draw$(style(i, j)) $(name(i)) -- $(name(j));"

diag(i::CartesianIndex{2}) = edge(i, i + CartesianIndex(1, 1))

antidiag(i::CartesianIndex{2}) = edge(i + CartesianIndex(1, 0), i + CartesianIndex(0, 1))

function grid(n::Int, m::Int, style=nostyle)
    nodes = [node(i) for i in CartesianIndices((n, m))]
    hedges = [edge(i, i + CartesianIndex(0, 1), style) for i in CartesianIndices((n, m-1))]
    vedges = [edge(i, i + CartesianIndex(1, 0), style) for i in CartesianIndices((n-1, m))]
    join(Iterators.flatten((nodes, hedges, vedges)), "\n")
end

function fig(fname::String, stuff...)
    open(fname, "w") do io
        println(io, "\\begin{tikzpicture}")
        for thing in stuff
            println(io, thing)
        end
        println(io, "\\end{tikzpicture}")
    end
end

fig(fname::String, diags::AbstractMatrix{Bool}) =
    fig(fname, grid((size(diags) .+ 1)...), [diag(i) for i in findall(diags)]...)

fig(fname::String, diags::AbstractMatrix{Bool}, antidiags::AbstractMatrix{Bool}) =
    fig(fname, grid((size(diags) .+ 1)...), [diag(i) for i in findall(diags)]...,
        [antidiag(i) for i in findall(antidiags)]...)

olabels(ribbons) =
    ["\\node at ($(i[2]-0.5), $(i[1]-0.5)) {$r};"
     for (r, indices) in enumerate(ribbons)
         for i in indices]

function onion(fname::String, ribbons)
    function style(i, j)
        if j[1] == i[1]
            i[1] > i[2] ? "[ultra thick]" : "[thin, gray]"
        else
            j[1] <= j[2] ? "[ultra thick]" : "[thin, gray]"
        end
    end
    n = length(ribbons) + 1
    fig(fname, grid(n, n, style), olabels(ribbons)...)
end

function figs()
    simple = eg.grow_simple(250+31)
    save("paper/simple_arcs.png", eg.arcs(
        simple, eg.Atom(1, 1), [250, 150, 100, 50, 25, 10]))
    save("paper/simple_ugly_arcs.png", eg.arcs(
        simple, eg.Atom(32, 16), [250, 150, 100, 50, 25, 10]))
    q8192 = deserialize("models/q8192")
    save("paper/nice_arcs.png", eg.arcs(
        q8192.diags, eg.Atom(32, 16), Ref(10) .* [250, 150, 100, 50, 25, 10], 10))
end

randiags(p, n, rng) = rand(rng, Float64, (n, n)) .< p

randarc(p, n, rng) = eg.sps(randiags(p,n,rng), true) .> n - 8

function randcirc(p, n, rng)
    g = eg.Grid(randiags(p, n, rng), randiags(p, n, rng))
    eg.sps(g, eg.Atom(div(n,2), div(n,2))) .>= div(n,2)-4
end

function randcircs()
    rng = eg.StableRNG(1)
    save("paper/images/randcircs3.png", hcat([randcirc(1/8+(i/15)/8, 512, rng) for i in 0:10]...))
end

using Images

using DataFrames, GLM

transp(bits) = [b ? RGBA(1, 1, 1, 0) : RGBA(0, 0, 0, 1) for b in bits]

function do_rmangle(rmangle)
    data = DataFrame(X=collect(0:255), Y=sum(.!rmangle, dims=2)[1:256])
    rmangle = transp(rmangle)
    for (x, y) in enumerate(predict(lm(@formula(Y ~ X), data), DataFrame(X=0:511)))
        rmangle[x, Int(round(y)):Int(round(y))+1] .= RGBA(1, 0, 0, 1)
    end
    rmangle
end


function sarc(x)
    n = size(x, 1) - 1
    sum(x[i] == (sqrt((i[1]-1)^2 + (i[2]-1)^2) > n - 8)
        for i in CartesianIndices(x))
end

function randarcs()
    rng = eg.StableRNG(1)
    n = 511
    save("paper/images/randarcs.png", transp(hcat([randarc(i/8, n, rng) for i in 0:8]...)))
    zdiags = [randiags(p, n, rng) for p in 0.13:0.01:0.21]
    zoomed = [eg.sps(z, true) .> n - 8 for z in zdiags]
    save("paper/images/randarcs_zoomed.png", transp(hcat(zoomed...)))
    @assert argmax(sarc.(zoomed)) == 5
    save("paper/images/best_randarc.png", transp(zoomed[5]))

    rmangle = eg.sps(zdiags[5], true)[:,1:128] .!= 0:511
    save("paper/images/best_randarc_mangle.png", transp(rmangle))
end

function simparcs()
    n = 511
    xoff = 32
    yoff = 64
    simple = eg.grow_simple(n + max(xoff, yoff))

    k = 128
    simple_leq = eg.grow_simple(k, true)
    save("paper/images/simple.png", colorview(Gray, [
        simple[i] == simple_leq[i] ? (1.0 - simple[i]) : 0.8 for i in eg.Atoms(k)]))

    ssps = eg.sps(simple[1:n,1:n], true)
    save("paper/images/best_simparc.png", transp(ssps .> n - 8))
    save("paper/images/best_simparc_mangle.png",
         transp(ssps[:,1:128] .!= 0:511))
    ssps = eg.sps(simple[1+xoff:n+xoff,1+yoff:n+yoff], true)
    save("paper/images/bad_simparc.png", transp(ssps .> n - 8))
    save("paper/images/bad_simparc_mangle.png",
         transp(ssps[:,1:128] .!= 0:511))
end

function comparcs(q512 = eg.grow_grid(512))
    n = 511
    ssps = eg.sps(q512.diags[1:n,1:n])

    x = BitMatrix(undef, n+1, n+1)
    for i in CartesianIndices(x)
        x[i] = sqrt((i[1]-1)^2+(i[2]-1)^2) > n - 8
    end

    save("paper/images/comparc.png", transp(ssps .> n - 8))
    save("paper/images/comparc_mangle.png",
         transp(ssps[:,1:128] .!= 0:511))

    z = sum((ssps .> n - 8) .!= x)

    ssps = eg.sps(q512.diags[32:n+31,64:n+63])
    save("paper/images/shift_comparc.png", transp(ssps .> n - 8))
    save("paper/images/shift_comparc_mangle.png",
         transp(ssps[:,1:128] .!= 0:511))

    z
end


function comparcs(diags::BitMatrix)
    n = 511
    ssps = eg.sps(diags[1:n,1:n])

    x = BitMatrix(undef, n+1, n+1)
    for i in CartesianIndices(x)
        x[i] = sqrt((i[1]-1)^2+(i[2]-1)^2) > n - 8
    end

    save("paper/images/comparc.png", transp(ssps .> n - 8))
    save("paper/images/comparc_mangle.png",
         transp(ssps[:,1:128] .!= 0:511))

    sum((ssps .> n - 8) .!= x)
end

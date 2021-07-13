module x

using CircularArrays: CircularArray

const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

TripleVec = Vector{Tuple{CartesianIndex{2}, Int}}

struct Torus
    n::Int
    k::Int
    d::Array{Int, 4}
    diags::CircularArray{Bool, 2, BitMatrix}
    triples::Matrix{TripleVec}
end

function triples_upto(n)::TripleVec
    ts = TripleVec()
    for (a, b, c) in [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25),
                      (20, 21, 29), (12, 35, 37), (9, 40, 41), (28, 45, 53),
                      (11, 60, 61), (16, 63, 65), (33, 56, 65), (48, 55, 73),
                      (13, 84, 85), (36, 77, 85), (39, 80, 89), (65, 72, 97)]
        b > n && continue
        a, b, c = [a, b, c] .* div(n, b)
        push!(ts, (CartesianIndex(a, b), c), (CartesianIndex(b, a), c))
    end
    ts
end

#for u in CartesianIndices(triples)
#    l = sqrt(sum(u.I.^2))
#    !isinteger(l) && continue


function Torus(n::Int, k::Int)
    d = repeat(CartesianIndices((k, k)) .|> x->sum(x.I), outer=(1, 1, n, n))
    triples = [TripleVec() for _ in onexy:(onexy * k)]
    for (u, l) in triples_upto(n)
        for j in onexy:u
            push!(triples[j], (u, l))
        end
    end
    Torus(n, k, d, CircularArray(falses(n, n)), triples)
end

Base.CartesianIndices(t::Torus) = CartesianIndices(t.diags)

wrap(t::Torus, i::CartesianIndex{2}) = CartesianIndex(mod.(i.I, axes(t.diags)))

function add_diag(t::Torus, u::CartesianIndex{2})::Nothing
    @assert !t.diags[u]
    t.diags[u] = true

    du = view(t.d, :, :, u)
    du[:, 1] = du[1, :] = 1:t.k
    diags = view(t.diags, u+onexy:u+(t.k - 1)*onexy)
    for i in CartesianIndices(diags)
        du[i+onexy] = 1 + (diags[i] ? du[i] : min(du[i+onex], du[i+oney]))
    end

    kxy = t.k * onexy
    for i in Iterators.take(CartesianIndices(du), t.k^2 - 1)
        v = wrap(t, u - kxy + i)
        dvu = maximum(i.I) == t.k ? t.k - minimum(i.I) : t.d[kxy - i, v]
        dv = view(t.d, onexy + kxy - i: kxy, v)
        dv[1] == dvu + 1 && continue
        dv .= min.(dv, dvu .+ view(du, onexy:i))
    end
    nothing
end

deltas(t::Torus) = [t.d[u, i] - l for i in CartesianIndices(t), (u, l) in t.triples[onexy]]

bound(t::Torus)::Int = maximum(abs.(deltas(t)))

function dvia(t::Torus, u::CartesianIndex{2}, v::CartesianIndex{2}, w::CartesianIndex{2})
    duv = minimum(v.I) == 1 ? maximum(v.I) : t.d[v - onexy, u] + 1
    dvw = minimum(w.I) == 0 ? maximum(w.I) : t.d[w, wrap(t, u + v)]
    duv + dvw
end

function addables(t::Torus, b::Int=bound(t))::Set{CartesianIndex{2}}
    ix = Set(i for i in CartesianIndices(t) if !t.diags[i])
    for i in CartesianIndices(t), j in onexy:(onexy * t.k)
        u = wrap(t, i + j - onexy)
        !(u in ix) && continue
        for (k, l) in t.triples[j]
            d_i_k = t.d[k, i] - 1
            d_i_k >= l - b  && continue
            if dvia(t, i, j, k - j) == d_i_k
                pop!(ix, u)
                break
            end
        end
    end
    ix
end

function score(t::Torus)
    b = bound(t)
    b = 0
    candidates = addables(t, b)
    scores = zeros(Int, t.n, t.n)
    #scores = rand(t.n, t.n)
    for i in CartesianIndices(t), j in onexy:(onexy * t.k)
        u = wrap(t, i + j - onexy)
        !(u in candidates) && continue
        for (k, l) in t.triples[j]
            d_i_k = t.d[k, i]
            dvia(t, i, j, k) >= d_i_k && continue
            scores[u] += d_i_k - l
        end
    end
    for a in candidates
        scores[a] += 1e6
    end
    scores
end

function tdist(n, u, v)
    x = abs.(v.I .- u.I)
    #sum(min.(x, n .- x))
    maximum(min.(x, n .- x))
end

function score(t::Torus, b)
    candidates = addables(t, b)
    scores = zeros(Int, t.n, t.n)
    for u in candidates
        c = deepcopy(t)
        add_diag(c, u)
        scores[u] = length(addables(c, b))
    end
    #for i in CartesianIndices(scores)
    #    scores[i] == 0 && continue
    #    scores[i] = scores[i] * 1e6 + tdist(t.n, i, last)
    #end
    scores
end

function fromint(n, i)
    t = Torus(n, n)
    for x in CartesianIndices(t)
        m = x[1] - 1 + (x[2] - 1) * n
        if (i >> m) % 2 == 1
            add_diag(t, x)
        end
    end
    t
end


brute(n) = (fromint(n, i) for i in 1:2^(n^2))

function add_all2(t::Torus, b)
    added = []
    while true
        scores = score(t, b)
        s = maximum(scores)
        s < 1 && return
        cs = [i for i in CartesianIndices(t) if scores[i] == s]
        for a in reverse(added)
            length(cs) == 1 && break
            tds = [tdist(t.n, c, a) for c in cs]
            mtd = maximum(tds)
            cs = [c for (c, td) in zip(cs,tds) if td == mtd]
        end
        ix = cs[1]
        println(s, " ", ix)
        add_diag(t, ix)
        push!(added, ix)
    end
    nothing
end

function add_all(t::Torus)
    for b in -1:1
        println("xxx ", b)
        #=
        while true
            s, ix = findmax(score(t, b))
            s < 1 && break
            println(s, " ", ix)
            add_diag(t, ix)
        end
        =#
        add_all2(t, b)
    end
    nothing
end

function exnil(n, k=n)
    t = Torus(n, n)
    for i in CartesianIndices(t)
        sum(i.I.%2) == 0 && add_diag(t, i)
    end
    add_all(t)
    t
end

function grow(n, z=4)
    t = Torus(z, z)
    add_all(t)
    for i in 1:n
        t2 = Torus(2*t.n, 2*t.n)
        for u in CartesianIndices(t)
            t.diags[u] && add_diag(t2, u*2)
        end
        t = t2
        add_all(t)
        println(i, " ", bound(t))
    end
    t
end



using Evolutionary
using Statistics

function frombits(xs)
    k = Int(sqrt(length(xs)))
    t = Torus(k, k)
    ix = CartesianIndices(t.diags)
    for i in eachindex(xs)
        if xs[i]
            add_diag(t, ix[i])
        end
    end
    t
end

function tobits(t)
    xs = falses(t.n^2)
    for i in 1:t.n^2
        xs[t.diags[i]


#maxabs(t) = maximum(abs.(err(t)))
#meanabs(t) = mean(abs.(err(t)))
#mse(t) = mean(err(t).^2)
#errstats(t) = (maxabs(t), meanabs(t), mse(t))

function oscore(xs)
    t = frombits(xs)
    se = ((t.d[u, i] - l)^2 for i in CartesianIndices(t), (u, l) in t.triples[onexy])

    1e6 * bound(t) + mean(se)
    #div(1e6 * maxabs(t), 1) * 1000 + meanabs(t)
    #mse(t)
end

function ga(k, inst=falses(k^2))
    Evolutionary.optimize(oscore, inst, GA(
        populationSize=1000,
        crossoverRate=0.5,
        mutationRate=0.02,
        epsilon=0.1,
        selection=susinv,
        crossover=discrete,
        mutation=flip), Evolutionary.Options(iterations=100, show_trace=true))
end






        for j in i:(i+u-onexy)
            t.diags[j] && continue
            if t.d[j - i + onexy, i]


        for (j, triples) in from(t, i)
            t.diags[j] && continue




        from

function foo()
    for i in unfill(t), (u, l) in t.triples
        t.d[u, i] <= l && continue

    for


function can_add_diag(t::Torus, u::CartesianIndex{2} b::Int=bound(t))::Bool
    t.diags[u] && return false
    for v in (u - onexy * t.k)
        t.d[

        for (w, l) in t.triples


    b = bound(t)
    [i for i in CartesianIndices(t) if feasible(t, i, b)]


function score(t::Torus, u::CartesianIndex{2})::Int
    dur = Vector{Int}(undef, t.m)
    dur[1] = dur[t.m] = t.k + 1
    dur[2:t.m-1] = 1 .+ view(t.d, view(t.lix, t.m-1:-1:2), wrap(t, u + onexy))

    ls = t.lix .|> l->wrap(t, u - l)
    dlu = Vector{Int}(undef, t.m)
    dlu[1] = dlu[t.m] = t.k
    dlu[2:t.m-1] = t.d[CartesianIndex.(view(t.lix, 2:t.m-1), view(ls, 2:t.m-1))]

    s = 0
    for j in eachindex(t.rix), i in eachindex(t.lix)
        diuj = dlu[i] + dur[j]
        dij = t.d[t.lix[i] + t.rix[j], ls[i]]
        if diuj < dij
            s += 3 * (diuj * dij - t.sed[i, j]) + 1
        end
    end
    s
end

function remove_diag(t::Torus, u::CartesianIndex{2})::Nothing
    @assert t.diags[u]
    t.diags[u] = false

    t2 = Torus(t.n, t.k)
    for v in CartesianIndices(t.diags)
        add_diag(t2, v)
    end
    copy!(t.dk, t2.dk)
    copy!(t.dm, t2.dm)
    nothing
end


function add_best(t::Torus, scores::AbstractMatrix)
    #d = density(t)
    #(s, _), u = findmax(collect(zip(scores, .-d)))
    s, u = findmax(scores)
    s <= 0 && return nothing
    if !t.diags[u]
        add_diag(t, u)
    else
        remove_diag(t, u)
    end
    for v in CartesianIndices(t.diags)
        scores[v] = score(t, v)
    end
    (s, u)
end

function add_all(t::Torus, scores::AbstractMatrix)
    for i in 1:length(t.diags)
        x = add_best(t, scores)
        #t = foo(t)
        x == nothing && return
        println(i, " ", x, " ", t.diags[x[2]])
    end
end

function err(t::Torus)
    de = CartesianIndices((t.m, t.m)) .|> i->sqrt(sum(i.I.^2))
    (t.d .- de)
end

function density(t::Torus)
    d = CircularArray(convert.(Float64, t.diags))
    for _ in 1:t.m
        for i in CartesianIndices(t.diags)
            #dlocal = (sum(d[j] for j in CartesianIndices(i-onexy:i+onexy)) - d[i]) / 8
            dlocal = (d[i - onexy] + d[i + onexy]) / 2
            d[i] = 0.9 * d[i] + 0.1 * dlocal
        end
    end
    d
end

function check(k)
    t = Torus(2*k+1, k)
    for i in CartesianIndices(t.diags)
        rand(Bool) && add_diag(t, i)
    end
    t, maximum(abs.(err(t)))
end



function expand(t::Torus)
    #big = Torus(2*t.n, div(2*t.n - 1, 4))
    big = Torus(t.w+4, t.k+1)
    for i in CartesianIndices(t.diags)
        if t.diags[i]
            add_diag(big, i)
            #add_diag(big, 2 * i)
        end
    end
    big
end

function grow(k, t = Torus(1, 0))
    scores = fill(score(t, onexy), t.n, t.n)
    add_all(t, scores)
    println(maximum(abs.(err(t))))
    #t = foo(t)
    for i in 2:k
        t = expand(t)
        println("n=",t.n)
        scores = [score(t, i) for i in CartesianIndices(t.diags)]
        add_all(t, scores)
        println(maximum(abs.(err(t))), " ", mean(abs.(err(t))))
        #t = foo(t)
    end
    t
end

using Statistics

function foo(t)
    z = mean(err(t).^2)
    t3 = t
    for i in CartesianIndices(t.diags)
        !t.diags[i] && continue
        t2 = Torus(t.n, t.k)
        for j in CartesianIndices(t.diags)
            t.diags[j] && j!=i && add_diag(t2, j)
        end
        if mean(err(t2).^2) < z
            t3 = t2
        end
    end
    t3
end

using Evolutionary
using Statistics

function frombits(xs)
    k = Int(div(sqrt(length(xs)), 2))
    t = Torus(2*k+1, k)
    ix = CartesianIndices(t.diags)
    for i in eachindex(xs)
        if xs[i]
            add_diag(t, ix[i])
        end
    end
    t
end

maxabs(t) = maximum(abs.(err(t)))
meanabs(t) = mean(abs.(err(t)))
mse(t) = mean(err(t).^2)

errstats(t) = (maxabs(t), meanabs(t), mse(t))

function oscore(xs)
    t = frombits(xs)
    #div(1e6 * maxabs(t), 1) * 1000 + meanabs(t)
    #mse(t)
    div(1e6 * maxabs(t), 1) * 1000 + mse(t)
end

function ga(k, inst=falses((2*k+1)^2))
    Evolutionary.optimize(oscore, inst, GA(
        populationSize=10000,
        crossoverRate=0.8,
        mutationRate=0.02,
        epsilon=0.1,
        selection=susinv,
        crossover=discrete,
        mutation=flip), Evolutionary.Options(iterations=10, show_trace=true))
end


#=

    ls = [wrap(t, u - i) for i in t.lix]

    dlu = Vector{Int}(undef, t.m)
    dlu[1] = dlu[t.m] = t.k
    dlu[2:t.m-1] = t.d[CartesianIndex.(view(t.lix, 2:t.m-1), view(ls, 2:t.m-1))]

    dut = Vector{Int}(undef, t.m)
    dut[1] = dut[t.m] = t.k
    dut[2:t.m-1] = t.d[view(t.tix, 2:t.m-1) .- [onexy], wrap(t, u + onexy)]

    s = 0
    for j in eachindex(t.tix), i in eachindex(t.lix)
        diuj = dlu[i] + dut[j] + 1
        dij = t.d[t.lix[i] + t.tix[j], ls[i]]
        if diuj < dij
            s += 3 * (diuj * dij - t.sed[i, j]) + 1
        end
    end
    s
end




function lindices(t::Torus, u::CartesianIndex{2})::Vector{CartesianIndex{4}}
    ixs = ((0:t.k-1,t.k:t.k), (t.k:t.k, t.k:-1:0))
    CartesianIndex.((i, wrap(t, u-i)) for ix in ixs for i in CartesianIndices(ix))
end

function tindices(t::Torus, u::CartesianIndex{2})::Vector{CartesianIndex{2}}
    ixs = ((k-1,2:i:

        0:t.k-1,t.k:t.k), (t.k:t.k, t.k:-1:0))
    CartesianIndex.((i, wrap(t, u-i)) for ix in ixs for i in CartesianIndices(ix))
end




                         i in Iterators.Flatten(

                        xindices(
    (i, k, u - CartesianIndex

    (0:k-1, k, u[1]-1:-1:u[1]-k+1, u[2]-k);
)

    dlu[2:k] = view(t.d, CartesianIndex.(1:k-1, k, u[1]-1:-1:u[1]-k+1, u[2]-k))
    dlu[k+1:t.m-1] = view(k, k:-1:1, u[1]-k,
    copy!
    copyto!(dlu, view(t.d, 1:k-1,
end


    mxy = uxy * t.m
    ll = o - mxy



        d_i_oll == d_

        d_i = view(t.d, :, :, wrapn(t, ll + i))
        d_i

        d_o_i = view(d_o, uxy:i)

        d_j = view(t.d, mxy + uxy - i:mxy, j)
        d_j_o = t.d[
        do_i



        j = foo
        du_i = view(t.d, uxy:i,

        d0 = minimum(i.I) == 1 ? (maximum(i.I) - 1) : d_o


        if i == uxy
            continue
        end
        du_o = view(d_o, uxy:i)
        du_i = view(t.d,

            view(d_o, i:t.m*uxy - i)





lowrapn(t::Torus, i::Int)::Int = i < 1 ? t.n - i : i
hiwrapn(t::Torus, i::Int)::Int = i > t.n ? i - t.n : i

lowrapn(t::Torus, i::Int, j::Int)::CartesianIndex{2} =
    CartesianIndex(lowrapn(t, i), lowrapn(t, j))

hiwrapn(t::Torus, i::Int, j::Int)::CartesianIndex{2} =
    CartesianIndex(hiwrapn(t, i), hiwrapn(t, j))

around(t::Torus, x::Int, y::Int) =
    ((i, hiwrapn(t, i + CartesianIndex(x, y)))
     for i in CartesianIndices((x:x+t.m-1, y:y+t.m-1)))




function remove_diag(t::Torus, o::CartesianIndex{2})
    for d_




function add_diag(t::Torus, x::Int, y::Int)
     _, dx, dy, dxy = CartesianIndices((0:1, 0:1))
    for (i, j) in around(t, x, y)
        d0 = minimum(i.I) == 1 ? (maximum(i.I) - 1) : t.d[i - dxy, j]
        d = view(t.d, i:CartesianIndex((t.k, t.k)), j)
        d_1 = view(d, : 1)
        copy!(d_1, min.(d_1, d0:d0+t.k-i[1]))
        d1_ = view(d, 1, :)
        copy!(d_1, min.(d_1, d0:d0+t.k-i[1]))
            d[1, :] = min.(d[1, :], d0:d0+j[2])
            for k in CartesianIndices(d)
                d[k + dxy] = 1 + (diags[k] ? d[k] : min(d[k+dx], d[k+dy]))
            end

        if
functio



sps_box(t::Torus, i::Int, j::Int, x::Int, y::Int) =
    view(t.dg, i:t.w, j:t.w, lowrapn(t, i + x - t.w, j + y - t.w))

sps_dgs(t::Torus, x::Int, y::Int) =
    ((sps_box(t, i, j, x, y), t.dg[i - 1, j - 1, x, y])
     for i in 2:t.w-1, j in 2:t.w-1)

function add_diag(t::Torus, x::Int, y::Int)
    for i in roi(t.n, x-t.w, x+t.w, y-t.w, y+t.w)
        for v in roi()
            for w in roi()



function update_dg(t::Torus, x::Int, y::Int)
    d_xy = view(t.dg, :, :, x, y)
    d_xy[1,:] = 1:t.n
    d_xy[:,1] = 1:t.n
    d_xy[2:t.w,2:t.w] = 1 .+ view(t.dg, 1:t.w-1, 1:t.w-1, hiwrapn(t, x+1, y+1))
    for (dur_ij, dxy_ij) in sps_dgs(t, x, y)
        dur_ij[1] == dxy_ij + 1 && continue
        copy!(dur_ij, min(dur_ij, dxy_ij .+ view(d_xy, size(dur_ij)...)))
    end
end

function dl(t::Torus, i::Int, x::Int, y::Int)
    (i == 1 || i == t.m) && return t.k
    i -= 1
    i <= t.k && return t.dg[i,t.k, x-i, y-t.k]
    i -= t.k
    t.dg[t.k,t.k-i,x-t.k,y-t.k+i]
end

function du(t::Torus, i::Int, x::Int, y::Int)
    (i == 1 || i == t.m) && return t.k
    i -= 1
    i <= t.k && return t.dg[t.k, i,x+1, y+1]
    i -= t.k
    t.dg[t.k-i,t.k,x+1,y+1]
end

function dlu(t::Torus, i::Int, j::Int, x::Int, y::Int)
    if i <= t.k
        dx, dy = i-1,t.k
    else
        dx, dy = t.k, t.m - i
    end
    x -= dx
    y -= dy
    if j <= t.k
        dx += t.k+1
        dy += j
    else
        dx += t.m - j
        dy += t.k
    end
    t.dg[dx, dy, x, y]
end



function score(t::Torus, x::Int, y::Int)
    s = 0
    for ,  in lx(t, x, y)
        if dvia() == t.dg[

            t.dg[x-lx,y-ly+1,lx,ly] + 1 + t.dg[ux-x+1,uy-y+1,ux,uy] <
            t.dg[

end



function update_scores(t::Torus, x::Int, y::Int)


function update_sps(t::Torus, x::Int, y::Int)
    at = view(t.sps, x, y)
    at[1,:] = 1:t.n
    at[:,1] = 1:t.n
    at[2:t.w,2:t.w] = 1 .+ view(t.sps, wrapn(t, x+1), wrapn(t, y+1),1:t.w-1,1:t.w-1)
    for (box, rx, ry) in boxes(t, x, y)

        for ur, d in ...
        ur = view(box, rx+1:end, ry+1:end)
        drxy = box[rx, ry]
        ur[1] == drxy + 1 && continue
        @. ur = min(ur, box[rx, ry] + view(at, shape(ur)...))
    end
end



        i in CartesianIndices((x:-1:x-t.w,y:y+t.w))
        ur = view(t.sps, i, x:end, y:end)
        ur .= min(ur, bb[dx, dy] +

    t.sps[x, y, 2:end, 2:end] = min(t.sps[x, y, 2:end, 2:end], 1 + t.sps[x - 1 1:end-1, y + 1:end-1]

    for (ll, ur, d) in llurd(t, x, y)
        ur .= min(ur, d



end


function flip(t::Torus, x::Int, y::Int)
    t.diags[x, y] = !t.diags[x, y]

end

const Score = Float64

function flip_best(t::Torus)::Union{Score, Nothing}
    s, i = findmax(t.scores)
    s <= 0 && return nothing
    flip(t, i...)
    s
end
=#

end
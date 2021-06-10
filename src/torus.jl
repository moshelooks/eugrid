const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

struct Torus
    n::Int
    k::Int
    m::Int
    w::Int
    ls::Vector{CartesianIndex{2}}
    sed::Matrix{Int}
    dk::Array{Int, 3}
    dm::Array{Int, 3}
    diags::CircularArray{Bool, 2, BitMatrix}
end

function Torus(n::Int, k::Int)
    m = 2*k + 1
    w = 2*m - 1
    @assert 1 <= w <= n
    ls = CartesianIndex.(0:m-1, m-1:-1:0)
    sed = ls .+ reshape(reverse(ls) .+ [onexy], 1, :) .|> x->sum(x.I.^2)
    dk = fill(2*k, (m, n, n))
    dm = fill(2*m, (w, n, n))
    Torus(n, k, m, w, ls, sed, dk, dm, CircularArray(falses(n, n)))
end

wrap(t::Torus, i::CartesianIndex{2}) = CartesianIndex(mod.(i.I, axes(t.diags)))

function score(t::Torus, u::CartesianIndex{2})
    if t.diags[u]
        return 0
        t2 = Torus(t.n, t.k)
        for v in Staircase(t, u)
            if t.diags[v] && v != u
                add_diag(t2, v)
            end
        end
        return -score(t2, u)
    end

    ls = t.ls .|> l->wrap(t, u - l)
    dlu = t.dk[CartesianIndex.(t.m:-1:1, ls)]
    dur = 1 .+ view(t.dk, :, wrap(t, u + onexy))

    s = 0
    for j in eachindex(dur), i in eachindex(dlu)
        diuj = dlu[i] + dur[j]
        dij = t.dm[t.m-i+j, ls[i]]
        if diuj < dij
            s += (3 * (diuj * dij - t.sed[i, j]) + 1)
            #s += 2 * (dij - sqrt(t.sed[i,j])) - 1
        end
    end
    #Int(floor(s*1000))
    s
end

struct Staircase
    t::Torus
    u::CartesianIndex{2}
end

function Base.iterate(I::Staircase, (i, j)=(0, 1-I.t.m))
    j == I.t.m && return nothing
    v = wrap(I.t, CartesianIndex(i, j) + I.u)
    i += 1
    if i + j == I.t.m
        i = 1 - I.t.m
        j += 1
    elseif i == I.t.m
        i = -I.t.m - j
        j += 1
    end
    (v, (i, j))
end

Base.eltype(::Type{Staircase}) = CartesianIndex{2}

Base.length(I::Staircase) = 3*I.t.m*(I.t.m - 1) + 1

function add_diag(t::Torus, u::CartesianIndex{2})::Nothing
    @assert !t.diags[u]
    t.diags[u] = true

    du = Matrix{Int}(undef, t.w, t.w)
    du[:, 1] = du[1, :] = 1:t.w
    diags = view(t.diags, u+onexy:u+(t.w - 1)*onexy)
    for i in UpperAntiTriangularIndices(t.w-2)
        du[i+onexy] = 1 + (diags[i] ? du[i] : min(du[i+onex], du[i+oney]))
    end
    t.dk[2:t.m-1, u] = antidiag(du, t.m - 2 - t.w)
    t.dm[:, u] = antidiag(du)

    du[:, t.w] = du[t.w, :] = t.w-1:-1:0
    diags = view(t.diags, u-(t.w-2)*onexy:u-onexy)
    for i in Iterators.reverse(LowerAntiTriangularIndices(t.w-2))
        j = i + onexy
        du[j] = 1 + (diags[i] ? du[j+onexy] : min(du[j+onex], du[j+oney]))
    end

    for i in Iterators.drop(UpperAntiTriangularIndices(t.w), 1)
        v = wrap(t, u+onexy-i)
        dvu = du[(t.w + 1)*onexy - i]

        dv = view(t.dm, i[2]:t.w-i[1]+1, v)
        dv .= min(dv, dvu .+ antidiag(du, 2 - sum(i.I)))

        dv = view(t.dk, i[2]:t.m-i[1]+1, v)
        dv .= min(dv, dvu .+ antidiag(du, 1 - t.k - sum(i.I)))
    end
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
    d = density(t)
    (s, _), u = findmax(collect(zip(scores, .-d)))
    #s, u = findmax(scores)
    s <= 0 && return nothing
    if !t.diags[u]
        add_diag(t, u)
    else
        remove_diag(t, u)
    end
    for v in Staircase(t, u)
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
    de = [sqrt(i^2+(t.w+1-i)^2) for i in 1:t.w]
    (t.dm .- de)
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

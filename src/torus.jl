
using CircularArrays

const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

struct Torus
    m::Int
    n::Int
    k::Int
    tix::Vector{CartesianIndex{2}}
    lix::Vector{CartesianIndex{2}}
    sed::Matrix{Int}
    d::Array{Int,4}
    diags::CircularArray{Bool, 2, BitMatrix}
end

function Torus(m::Int, n::Int)
    @assert 1 <= m <= n
    @assert m % 2 == 1
    k = div(m, 2)
    ixs = ((k+1:k+1,1:k), (k+1:-1:1,k+1:k+1))
    tix = CartesianIndex.(i for ix in ixs for i in CartesianIndices(ix))
    lix = CartesianIndex.(i - onexy for i in reverse(tix))
    sed = lix .+ reshape(tix, 1, :) .|> x->sum(x.I.^2)
    d = repeat(CartesianIndices((m, m)) .|> x->sum(x.I), outer=(1, 1, n, n))
    Torus(m, n, k, tix, lix, sed, d, CircularArray(falses(n, n)))
end

wrapn(t::Torus, i::CartesianIndex{2})::CartesianIndex{2} =
    CartesianIndex(mod.(i.I, axes(t.diags)))

function add_diag(t::Torus, u::CartesianIndex{2})::Nothing
    @assert !t.diags[u]
    t.diags[u] = true

    du = view(t.d, :, :, u)
    du[:, 1] = du[1, :] = 1:t.m
    diags = view(t.diags, u+onexy:u+(t.m - 1)*onexy)
    for i in CartesianIndices(diags)
        du[i+onexy] = diags[i] ? du[i] : min(du[i+onex], du[i+oney])
    end

    mxy = t.m * onexy
    for i in Iterators.take(CartesianIndices(du), t.m^2 - 1)
        v = wrapn(t, u - mxy + i)
        dv = view(t.d, :, :, v)
        dvu = maximum(i.I) == t.m ? t.m - minimum(i.I) : dv[mxy - i]

        dv = view(dv, onexy + mxy - i: mxy)
        dv[1] == dvu + 1 && continue
        #println(size(dv))
        #println(size(dvu))
        #println(size(view(du, onexy:i)))
        dv .= min.(dv, dvu .+ view(du, onexy:i))
    end
    nothing
end

function score(t::Torus, u::CartesianIndex{2})::Int
    ls = [wrapn(t, u - i) for i in t.lix]

    dlu = Vector{Int}(undef, t.m)
    dlu[1] = dlu[t.m] = t.k
    dlu[2:t.m-1] = t.d[CartesianIndex.(view(t.lix, 2:t.m-1), view(ls, 2:t.m-1))]

    dut = Vector{Int}(undef, t.m)
    dut[1] = dut[t.m] = t.k
    dut[2:t.m-1] = t.d[view(t.tix, 2:t.m-1) .- [onexy], wrapn(t, u + onexy)]

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



#=
function lindices(t::Torus, u::CartesianIndex{2})::Vector{CartesianIndex{4}}
    ixs = ((0:t.k-1,t.k:t.k), (t.k:t.k, t.k:-1:0))
    CartesianIndex.((i, wrapn(t, u-i)) for ix in ixs for i in CartesianIndices(ix))
end

function tindices(t::Torus, u::CartesianIndex{2})::Vector{CartesianIndex{2}}
    ixs = ((k-1,2:i:

        0:t.k-1,t.k:t.k), (t.k:t.k, t.k:-1:0))
    CartesianIndex.((i, wrapn(t, u-i)) for ix in ixs for i in CartesianIndices(ix))
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

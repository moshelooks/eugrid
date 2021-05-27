
struct Torus
    diags::BitMatrix
    scores::Matrix{Float64}
    dg::Array{Int16,4}
    de::Matrix{Float64}
end

lowrapn(t::Torus, i::Int)::Int = i < 1 ? t.n - i : i
hiwrapn(t::Torus, i::Int)::Int = i > t.n ? i - t.n : i

lowrapn(t::Torus, i::Int, j::Int)::CartesianIndex{2} =
    CartesianIndex(lowrapn(t, i), lowrapn(t, j))

hiwrapn(t::Torus, i::Int, j::Int)::CartesianIndex{2} =
    CartesianIndex(hiwrapn(t, i), hiwrapn(t, j))

sps_box(t::Torus, i::Int, j::Int, x::Int, y::Int) =
    view(t.dg, i:t.w, j:t.w, lowrapn(t, i + x - t.w, j + y - t.w))

sps_dgs(t::Torus, x::Int, y::Int) =
    ((sps_box(t, i, j, x, y), t.dg[i - 1, j - 1, x, y])
     for i in 2:t.w-1, j in 2:t.w-1)

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

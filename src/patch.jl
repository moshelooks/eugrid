
function shortest_paths(diags::AbstractMatrix{Bool})::Matrix{Int}
    n, m = size(diags)
    d = Matrix{Int}(undef, n+1, m+1)
    d[1:n+1,1] = 1:n+1
    d[1, 1:m+1] = 1:m+1
    _, dx, dy, dxy = CartesianIndices((0:1, 0:1))
    for i in CartesianIndices(diags)
        d[i + dxy] = 1 + diags[i] ? d[i] : min(d[i+dx], d[j+dy])
    end
    d
end

const IntD = Int16

struct Patch
    dil::Array{IntD,3}
    dit::Array{IntD,3}
end

function graph_distances(p::Patch)
    n = size(p.dil, 1)
    k = div(n, 2)
    dg = Matrix{IntD}(undef, n, n)
    copyto!(dg, CartesianIndices((k,n)), p.dit, CartesianIndices((k+1:n-1,1,n:n)))
    copyto!(dg, CartesianIndices((k+1:n,n)), p.dit, CartesianIndices((n:n, 1:k+1, n)))
    dg
end


function score()
    s = 0

#    for i,
    #        if di + dj + 1 <
end


function add_diag(t::Torus, p::Patch, x::Int, y::Int)
    tix = max(1,y-t.k):min(t.n,x+t.k)
    dxyt = g.dit[x,y+1,tix]
    sps = shortest(view
    dit = view(p.dit, x+1:n+1, y:-1:1, min(1,y-k):min(n,x+k))

    for i in CartesianIndices((x+1:n+1,1:y))



function score_add(diags, x, y)
    l = sps(view(diags))
    u = sps(view(diags))





struct Torus

end

invalidated(t::Torus, i::CartesianIndex{2}) =
    wrap(t).CartesianIndices(i




function Patch(n::Int)::Patch


function foo()


function remove_diag(p::Patch, x::Int, y::Int)

function add_diag(p::Patch, x::Int, y::Int)
    lix = x:n+y-1
    dxyl = g.dil[x, y, lix]
    if x > 1 && y < n

    p.dil[x:end,y:end,

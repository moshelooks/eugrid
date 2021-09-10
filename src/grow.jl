const Atom = CartesianIndex{2}
const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

function choose_diag!(a::Atom, tl, l, t, br)::Bool
    n = length(l)
    br[1] = l[1] + 1
    br[n+1] = t[n] + 1
    br = view(br, 2:n)
    br .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1

    s = 0.0
    for i in 1:a[1]
        de = sqrt(i^2 + a[2]^2)
        s += abs(de - br[i])
        s -= abs(de - tl[i] - 1)
    end
    for i in a[1]+1:a[1]+a[2]-1
        de = sqrt(a[1]^2 + (a[1]+a[2]-i)^2)
        s += abs(de - br[i])
        s -= abs(de - tl[i] - 1)
    end

    if s > 0
        br .= view(tl, 1:n-1) .+ 1
        true
    else
        false
    end
end

function choose_children!(diags::BitMatrix, a::Atom, grandparents, parents, children)
    Threads.@threads for j in 1:size(children)[2]
        tl = view(grandparents, :, j)
        l = view(parents, :, j)
        t = view(parents, :, j+1)
        br = view(children, :, j)
        a_j = a + CartesianIndex(-1, 1) * (j - 1)
        diags[a_j] = choose_diag!(a_j, tl, l, t, br)
    end
end

function grow_diags(n::Int)::BitMatrix
    diags = BitMatrix(undef, n, n)
    buffer = Array{Int, 3}(undef, 2*n+1, n+2, 3)
    grandparents = view(buffer, 1:1, 1:1, 1)
    fill!(grandparents, 0)
    parents = view(buffer, 1:2, 1:2, 2)
    parents[1, 1] = parents[2, 2] = 0
    parents[1, 2] = parents[2, 1] = 1
    children = view(buffer, 1:3, 1:3, 3)
    for i in 1:n
        children = view(buffer, 1:i+2, 1:i+2, mod1(i+2, 3))
        children[:, 1] .= 0:i+1
        choose_children!(diags, Atom(i, 1), grandparents, parents, view(children, :, 2:i+1))
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)

    for i in n-1:-1:1
        children = view(buffer, 1:2*n-i+2, 1:i, mod1(2*n-i+2, 3))
        choose_children!(diags, Atom(n, n-i+1), grandparents, parents, children)
        grandparents = view(parents, :, 2:size(parents)[2])
        parents = children
    end
    diags
end

function score(g)
    g = BitMatrix(g)
    n = size(g)[1]
    step = Int(n / 8)
    sum(eg.arc(g, step) .!= eg.earc(n, step)) / n^2
end

function grow_simple(n::Int)::BitMatrix
    ds = Matrix{Int}(undef, n+1, n+1)
    ds[:, 1] = ds[1, :] = 0:n
    diags = BitMatrix(undef, n, n)
    for i in CartesianIndices(diags)
        tl = ds[i]
        l = ds[i+onex]
        t = ds[i+oney]
        de = sqrt(i[1]^2 + i[2]^2)
        d_diag = tl + 1
        d_nodiag = min(l, t) + 1
        if abs(de - d_diag) < abs(de - d_nodiag)
            ds[i + onexy] = d_diag
            diags[i] = true
        else
            ds[i + onexy] = d_nodiag
            diags[i] = false
        end
    end
    diags
end

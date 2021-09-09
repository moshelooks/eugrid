const Atom = CartesianIndex{2}
const onex, oney, onexy = CartesianIndices((0:1, 0:1))[2:end]

function grow2(tl, l, t, a, bcell)
    n = length(l)
    bcell[1] = l[1] + 1
    bcell[n+1] = t[n] + 1
    bcell = view(bcell, 2:n)
    bcell .= min.(view(l, 2:n), view(t, 1:n-1)) .+ 1

    s = 0.0
    for i in 1:a[1]
        de = sqrt(i^2 + a[2]^2)
        s += abs(de - bcell[i])
        s -= abs(de - tl[i] - 1)
    end
    for i in a[1]+1:a[1]+a[2]-1
        de = sqrt(a[1]^2 + (a[1]+a[2]-i)^2)
        s += abs(de - bcell[i])
        s -= abs(de - tl[i] - 1)
    end

    if s > 0
        bcell .= view(tl, 1:n-1) .+ 1
        true
    else
        false
    end
end

function grow2(n)
    diags = BitMatrix(undef, n, n)
    grandparents = zeros(Int, 1, 1)
    parents = [0 1; 1 0]
    for i in 1:n
        children = Matrix{Int}(undef, i+2, i+2)
        children[:, 1] .= 0:i+1
        a = Atom(i, 1)
        for j in 2:i+1
            tl = view(grandparents, :, j-1)
            l = view(parents, :, j-1)
            t = view(parents, :, j)
            diags[a] = grow2(tl, l, t, a, view(children, :, j))
            a += CartesianIndex(-1, 1)
        end
        children[:, i+2] .= i+1:-1:0
        grandparents = parents
        parents = children
    end
    grandparents = view(grandparents, :, 2:n)
    parents = view(parents, :, 2:n+1)
    for i in n-1:-1:1
        children = Matrix{Int}(undef, 2*n-i+2, i)
        a = Atom(n, n-i+1)
        for j in 1:i
            tl = view(grandparents, :, j)
            l = view(parents, :, j)
            t = view(parents, :, j+1)
            diags[a] = grow2(tl, l, t, a, view(children, :, j))
            a += CartesianIndex(-1, 1)
        end
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

function addable(t::Torus, b::Int, u::CartesianIndex{2})::Bool
    t.diags[u] && return false
    for j in onexy:(onexy * t.k)
        i = wrap(t, u - j + onexy)
        for (k, l) in t.triples[j]
            d_i_k = t.d[k, i] - 1
            d_i_k >= l - b  && continue
            dvia(t, i, j, k - j) == d_i_k && return false
        end
    end
    true
end

struct Score
    t::Torus
    u::CartesianIndex{2}
    addables::Vector{CartesianIndex{2}}
end

Score(t::Torus, b::Int)::Score =
    Score(t, CartesianIndex(0, 0), filter(u->addable(t,b,u), CartesianIndices(t)))

using Statistics

function serr(t)
    se = ((t.d[u, i] - l)^2 for i in CartesianIndices(t), (u, l) in t.triples[onexy])
    #se = (abs(t.d[u, i] - l) for i in CartesianIndices(t), (u, l) in t.triples[onexy])
    mean(se)
end

function score(t::Torus, u::CartesianIndex)
    t = deepcopy(t)
    add_diag(t, u)
    serr(t)
end

function add_all2(t::Torus, b::Int)::Torus
    addables = filter(u->addable(t,b,u), CartesianIndices(t))
    s0 = serr(t)
    while !isempty(addables)
        bound(t) <= b && break
        s, i = findmin([score(t, u) for u in addables])
        s >= s0 && break
        s0 = s
        u = addables[i]
        add_diag(t, u)
        addables = filter(u->addable(t,b,u), CartesianIndices(t))
        println(s, " ", length(addables), " ", u, " ", bound(t))
    end
    t
end


#=
function score(prev::Score, b::Int, u::CartesianIndex, best::Int)
    t = deepcopy(prev.t)
    add_diag(t, u)
    #=
    t = prev.t
    d = Array{Int, 4}(undef, t.k, t.k, t.k, t.k)
    kxy = t.k * onexy
    du = view(d, :, :, kxy)
    du[:, 1] = du[1, :] = 1:t.k
    diags = view(t.diags, u+onexy:u+(t.k - 1)*onexy)
    for i in CartesianIndices(diags)
        du[i+onexy] = 1 + (diags[i] ? du[i] : min(du[i+onex], du[i+oney]))
    end

    for i in Iterators.take(CartesianIndices(du), t.k^2 - 1)
        v = wrap(t, u - kxy + i)
        dvu = maximum(i.I) == t.k ? t.k - minimum(i.I) : t.d[kxy - i, v]
        dv = view(d, onexy + kxy - i: kxy, i)
        dv[1] == dvu + 1 && continue
        dv .= min.(dv, dvu .+ view(du, onexy:i))
    end
    =#
    addables = Vector{CartesianIndex{2}}()
    for (n, i) in enumerate(prev.addables)
        addable(t, b, i) && push!(addables, i)
        length(addables) + length(prev.addables) - n <= best + 1 && return nothing, best
    end
    Score(t, u, addables), length(addables)
end

function add_all2(t::Torus, b::Int)::Torus
    state = Score(t, b)
    while true
        next, best = state, 0
        symmetries = Set{BitMatrix}()
        for u in state.addables
            sym = symmetrize(state.t, u)
            sym in symmetries && continue
            transpose(sym) in symmetries && continue
            push!(symmetries, sym)
            next_u, best = score(state, b, u, best)
            if !isnothing(next_u)
                next = next_u
                best == length(state.addables) -1 && break
            end
        end
        best < 1 && break
        state = next
        print(best, " ", state.u, " ")
        println(bound(state.t))
    end
    state.t
end
=#

function add_all(t::Torus, bmax::Int=1)::Torus
    for b in -1:bmax
        println("xxx ", b, " ", extrema(deltas(t)))
        t = add_all2(t, b)
    end
    t
end

function modbound(n::Int)::Int
    maximum(a/(a+b-c) for (a,b,c) in Pythagorean.scaled_triples(n))
end

function mind(x::Int, y::Int, m::Int)::Int
    x, y = minmax(x, y)
    max(y, x + div((m-1)*y-x+1, m))
end

function tribound(n::Int)::Int
    n in 10:11 && return 3
    m = 1
    for (a, b, c) in Pythagorean.scaled_triples(n)[1:2:end]
        while mind(a, b, m) < c
            m += 1
        end
    end
    m
end

function init(n::Int, b, m::Int=tribound(n))::Torus
    t = Torus(n)
    for i in 1:b:(n-b+1)
        for j in 1:n
            #u = wrap(t, CartesianIndex(n-j+1,i+j-1))
            u = CartesianIndex(i, j)
            add_diag(t, u)
        end
    end
    for i in 1:m:(n-m+1)
        for j in 1:n
            #u = wrap(t, CartesianIndex(n-j+1,i+j-1))
            u = CartesianIndex(j, i)
            t.diags[u] && continue
            add_diag(t, u)
        end
    end
    t
end

function init2(n::Int, m, k, kk, ll)::Torus
    t = init(n, m)
    for i in (k+1):m:(n-m+1+k)
        for j in 1:ll:kk
            u = wrap(t, CartesianIndex(n-j+1,i+j-1))
            add_diag(t, u)
        end
    end
    t
end

function init3(n::Int, m, kk, k, z)::Torus
    t = init2(n, m, kk)
    for i in (k+1):m:(n-m+1+k)
        for j in 1:z:(n-m+1+z)
            u = wrap(t, CartesianIndex(n-j+1,i+j-1))
            add_diag(t, u)
        end
    end
    t
end


function check_mind(m::Int)
    t = init(m * (m+1), m)
    for i in CartesianIndices((m, m))
        @assert minimum(t.d[i, :, :]) == mind(i.I..., m)
    end
end



function exnil(n::Int, k::Int=n, bmax::Int=1)::Torus
    #t = Torus(n, k)
    #=
    for i in CartesianIndices(t)
        #sum(i.I.%2) == 0 && add_diag(t, i)
        #i[1]%4==0 && i[2]%4!=0 && add_diag(t, i)
        i[1] + i[2] == t.n + 1 && add_diag(t, i)
    end
    =#
    #=
    for i in 1:n
        for j in 1:n
            t2 = deepcopy(t)
            u = wrap(t, CartesianIndex(n-j+1,i+j-1))
            add_diag(t2, u)
            minimum(deltas(t2)) >= 0 && (t = t2)
        end
    end
    =#
    #=
    for u in CartesianIndices(t)
        t2 = deepcopy(t)
        add_diag(t2, u)
        minimum(deltas(t2)) >= 0 && (t = t2)
    end
    =#
    t = init(n)

    add_all(t, bmax)
end

function gx(d::Int)::Torus
    t = exnil(4, 4, -1)
    for i in 1:d
        t2 = Torus(2*t.n, 2*t.n)
        for u in CartesianIndices(t)
            if t.diags[u] && maximum(u.I) <= t.n - i + 1 && minimum(u.I) > i - 1
                #add_diag(t2, u*2-onexy)
                #add_diag(t2, u*2-oney)
                add_diag(t2, u)
                add_diag(t2, u + onex * (t.n))
                add_diag(t2, u + oney * (t.n))
                add_diag(t2, u + onexy * (t.n))
            end
            #=
            v = CartesianIndex(u[2],u[1])
            if t.diags[v]
                add_diag(t2, v + onex * (t.n))
                add_diag(t2, v + oney * (t.n))
            end
            =#
        end
        t = t2
        println("yyy ", i, " ", extrema(reshape(deltas(t), :)))
        #t = add_all(t, i == d ? 1 : max(-1,i-d+1))
        t = add_all(t, i == d ? 1 : 0)
    end
    t
end


function grow(k::Int,l::Int=1)
    for n in 4:k
        t = exnil(n,k)
        b = bound(t)
        println()
        println("NB ", n," ",b)
        b <= l && return
    end
    nothing
end

function circ(t::Torus, o::CartesianIndex{2}=CartesianIndex(1, 1), r::Int=t.k-1)
    x = t.d[:, :, o] .<= r
    [rot180(x) rotl90(x); rotr90(x) x]
end



#=

struct Adder
    t::Torus
    b::Int
    addables::Vector{CartesianIndex{2}}
end

Adder(t::Torus, b::Int)::Adder = Adder(t, b, filter(u->addable(t,b,u), CartesianIndices(t)))

function add_diag(t::Torus, u::CartesianIndex{2}, a::Adder)::Nothing
    add_diag(t, u)
    filter!(u->addable(t,b,u), a.addables)
    nothing
end

function score(a::Addable)::Matrix{Int}
    scores = zeros(Int, a.t.n, a.t.n)

    kxy = t.k * onexy
    for u in a.addables
        du = t.d[ :, :, u]
        du[:, 1] = du[1, :] = 1:t.k
        diags = view(t.diags, u+onexy:u+(t.k - 1)*onexy)
        for i in CartesianIndices(diags)
            du[i+onexy] = 1 + (diags[i] ? du[i] : min(du[i+onex], du[i+oney]))
        end

        dvs = [du]

        for i in Iterators.take(CartesianIndices(du), t.k^2 - 1)
            v = wrap(t, u - kxy + i)
            dvu = maximum(i.I) == t.k ? t.k - minimum(i.I) : t.d[kxy - i, v]
            dv = t.d[onexy + kxy - i: kxy, v]
            dv[1] == dvu + 1 && continue
            dv .= min.(dv, dvu .+ view(du, onexy:i))
            push!(dvs, dv)
        end

        for i in CartesianIndices(du)
            v = wrap(t, u - kxy + i)



        c = deepcopy(t)
        add_diag(c, u)
        scores[u] = count(addable(c, b, i) for i in a.addables)
    end
    scores
end
=#


function fromint(n, x)
    t = Torus(n, n)
    for m in 0:(n-1)
        if (x >> m) % 2 == 1
            i = m+1
            for j in 1:n
                u = wrap(t, CartesianIndex(n-j+1,i+j-1))
                #=t2 = deepcopy(t)
                add_diag(t2, u)
                minimum(deltas(t2)) >= -1 && (t = t2)=#
                add_diag(t, u)
            end
        end
    end
    t
end


brute(n) = (fromint(n, i-1) for i in 1:2^n)

function blog(n)
    m = n^2
    for (i, t) in enumerate(brute(n))
        b = bound(t)
        if b < m
            println(b)
            m = b
        end
        if i%1000 == 0
            println(i, " (", div(i*100, 2^n), "%)")
        end
    end
    m
end

function rlog(n)
    m = n^2
    i = 1
    while m > 0
        x = rand(1:2^n)-1
        t = fromint(n, x)
        a, b = extrema(deltas(t))
        if a >= -1 && b < m
            println(i, " ", (a, b), " ", x)
            m = b
        end
        i+=1
        if i%1000 == 0
            println(i)
        end
    end
    m
end

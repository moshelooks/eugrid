function sps(diags::AbstractMatrix{Bool})::Matrix{Int}
    d = Matrix{Int}(undef, size(diags) .+ 1)
    d[:, 1] = 0:size(diags)[1]
    d[1, :] = 0:size(diags)[2]
    for i in CartesianIndices(diags)
        d[i + onexy] = 1 + (diags[i] ? d[i] : min(d[i+onex], d[i+oney]))
    end
    d[2:end, 2:end]
end

function iseugrid(diags::AbstractMatrix{Bool},
                  ts::Vector{Triple}=all_triples(maximum(size(diags))))::Bool
    for i in CartesianIndices(diags)
        d = sps(diags[i:CartesianIndex(size(diags))])
        x, y = size(d)
        for t in ts
            if t.a <= x && t.b <= y
                d[t.a, t.b] != t.c && return false
            end
        end
    end
    true
end

function fromint(n, x)::BitMatrix
    diags = BitMatrix(undef, n, n)
    for m in 0:n^2-1
        diags[m+1] = (x >> m) % 2 == 1
    end
    diags
end


brute(n) = (fromint(n, i-1) for i in Int128(1):Int128(2)^n^2)

function trace(kernel, basis, dmat, dbg=false)
    dset = Set(findall(dmat))
    o = Onion(kernel, basis, only(Set(size(dmat))))
    for ribbon in o.ribbons
        diags = intersect(dset, ribbon)
        dbg && println(diags)
        while true
            diags == popfirst!(o.solvers[end]) && break
        end
        didwrap = wrap!(o, diags)
        @assert didwrap
        if dbg && !iscomplete(o)
            a = only(o.solvers[end].stack)
            println(a.clauses)
            println(a.free)
            println(a.affirmed)
            println()
        end
    end
end

function trace_blocked(kernel, basis, dmat, dbg=false)
    dset = Set(findall(dmat))
    o = Onion(kernel, basis, only(Set(size(dmat))))
    for ribbon in o.ribbons
        diags = intersect(dset, ribbon)
        dbg && println(diags)
        while true
            #isempty(o.solvers[end]) && return o
            seen = popfirst!(o.solvers[end])
            diags == seen && break
            if length(o.solvers) < length(o.ribbons)
                ds = exterior_distances(o.membranes[end], seen)
                cs = constraints(o.box, ds)
                #println(cs)
                if !issatisfiable(cs)
                    #println("XXX ", violations(cs))
                    #=for v in violations(cs)
                        if v.u == Atom(4, 3)
                            return cs
                        end
                    end=#
                    bs = blockers(cs, o.membranes[end].interior, seen)
                    block!(o.solvers[end], bs)
                    for (b, v) in zip(bs, violations(cs))
                        if isblocked(Assignment([], Set(), diags), b)
                            println(b)
                            println(v)
                            println(seen)
                            return cs
                        end
                    end
                end
            end
        end
        didwrap = wrap!(o, diags)
        @assert didwrap
        if dbg && !iscomplete(o)
            a = only(o.solvers[end].stack)
            #println(a.clauses)
            #println(a.free)
            #println(a.affirmed)
            println()
        end
    end
end

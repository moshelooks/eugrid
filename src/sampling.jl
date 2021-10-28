"""
The basic idea here is to separate out the generation of random variates from the
computational heavy lifting so that the latter can be multi-threaded and we'll still have
reproducible results regardlessof how many threads there are.

10,000 xs for eccentricity

~ 0.35 * n * 100 ys for geocount

"""

struct Sample
    seed::Int

    x::Vertex
    y::Vertex
    z::Vertex

    eccentricity_x::Int
    euclidean_eccentricity_x::Float64

    dxy::Int
    geocount_xy::Int

    c::Vertex
    dcz::Int
end

function Sample(g::Grid, seed::Int, dxy::Int)::Sample
    rng = StableRNG(seed)

    x = rand(rng, vertices(g)[onexy*(dxy+1):onexy*(g.n-dxy)])
    dx = sps(g, x)

    eccentricity_x = eccentricity(g, x, dx)
    euclidean_eccentricity_x = euclidean_eccentricity(g.n, x)

    dxy == 0 && return Sample(
        seed, x, x, x, eccentricity_x, euclidean_eccentricity_x, dxy, 1, x, 0)

    y = rand(rng, circle_points(g, x, dxy, dx))
    @assert y != x

    dy = sps(g, y, dxy + 1)
    geocount_xy = sum(geodesics(g, x, y, dx, dy))

    dxy < (g.n - 1) / 16 && return Sample(
        seed, x, y, y, eccentricity_x, euclidean_eccentricity_x, dxy, geocount_xy, y, 0)

    z = rand(rng, two_circle_points(g, x, y, dxy, dx, dy, strict=false))
    @assert z != y

    c = rand(rng, midpoints(g, x, y, dx, dy))
    dcz = d2i(distance(g, c, z), RoundDown)

    Sample(seed, x, y, z, eccentricity_x, euclidean_eccentricity_x, dxy, geocount_xy, c, dcz)
end

delta_eccentricity(samples)::Vector{Float64} =
    [s.eccentricity_x - s.euclidean_eccentricity_x for s in samples]

function geodesic_exponent(samples)::Float64
    samples = [s for s in samples if s.x != s.y]
    xs = [log2(s.dxy) for s in samples]
    ys = [log2(s.geocount_xy) for s in samples]
    GLM.coef(GLM.lm(@GLM.formula(Y ~ X), DataFrames.DataFrame(X=xs, Y=ys)))[2]
end

delta_sqrt3(samples)::Vector{Float64} =
    [abs(sqrt(3) - 2s.dcz / s.dxy) for s in samples if s.y != s.z]

function score(samples)
    de = delta_eccentricity(samples)
    mu_de = Statistics.mean(de)
    ds = delta_sqrt3(samples)
    [mu_de, mu_de + Statistics.std(de),
     geodesic_exponent(samples),
     isempty(ds) ? NaN : Statistics.mean(ds)]
end

function sample(g::Grid, k::Int; seed::Int=1, min_geo::Int=4, max_geo::Int=2*div(7*g.n, 40))
    @assert max_geo >= min_geo
    nzeros = max(0, round(Int, k * (100 / (max_geo - min_geo + 1) - 1), RoundUp))
    dxys = repeat([zeros(Int, nzeros); min_geo:max_geo], k)
    samples = Vector{Sample}(undef, length(dxys))
    Threads.@threads for i in eachindex(dxys, samples)
        samples[i] = Sample(g, i + (seed - 1) * length(dxys), dxys[i])
    end
    samples
end

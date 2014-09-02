module CMSketch

type Count
    # sketch parameters
    epsilon::Float64
    delta::Float64
    # parameters for count array
    width::Int
    depth::Int
    # array of counts
    count::Array{Int,2}
    # dimension of vector
    N::Int
    # prime number
    P::Int
    # vector of (a, b)
    hash_params::Vector{(Int, Int)}

    function Count(epsilon, delta, N, P)
        width = iceil(exp(epsilon))
        depth = iceil(log(1/delta))
        count = zeros(depth, width)

        # N is the dimension of the input data vector
        N <= P || error("N < P!")
        # P must be prime
        isprime(P) || error("P needs to be prime.")

        # generate parameters for D=depth hash functions
        a = rand(1:P-1, depth)
        b = rand(0:P-1, depth)
        hash_params = collect(zip(a, b))

        new(epsilon, delta, width, depth, count, N, P, hash_params)
    end
end

# 2-universal hash family (adjusted for 1-indexing)
function hash_fn(a::Int, b::Int, width::Int, P::Int, x::Int)
    x -= 1
    mod(mod(a*x+b, P), width) + 1
end

# update the count matrix given (it, ct) stream
function update(sketch::Count, it::Int, ct::Int)

    for (j, (a,b)) in enumerate(sketch.hash_params)
        ind = hash_fn(a, b, sketch.width, sketch.P, it)
        sketch.count[j, ind] += ct
    end
end

# approximate a point query Q(i)
function point_query(sketch::Count, i::Int)

    function get_counts(j)
        (a,b) = sketch.hash_params[j]
        ind = hash_fn(a, b, sketch.width, sketch.P, i)
        sketch.count[j, ind]
    end

    counts = [get_counts(j) for j in 1:sketch.depth]
    minimum(counts)
end

end

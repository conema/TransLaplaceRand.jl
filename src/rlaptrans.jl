function rlaptrans(n::Integer, ltpdf::Function, maxiters::Integer=500, x0::Number=1.0, xinc::Number=2.0)
    # Generate random uniform number for inversion method
    U = sort(rand(Uniform(0, 1), n))

    # Definition of inverted pdf and cdf
    cdf = Fourier_series(t -> ltpdf(t)/t)
    pdf = Fourier_series(t -> ltpdf(t))

    # Calculate first pdf and cdf to be use as starting value for root algorithm
    t0 = x0

    # Find upper bound for root algorithm 
    ti = t0

    ctr = 0
    while ctr <= maxiters && cdf(ti) < U[n]
        ti *= xinc

        ctr += 1
    end

    ctr >= maxiters && error("Upper bound not found")

    # Root algorithm
    lb = 0.0
    ub = ti

    pdf2(t, u) = pdf(t)
    cdf2(t, u) = cdf(t)-u

    xrand = Vector{Float64}(undef, n)

    for (i, u) in enumerate(U)
        ran = find_zero(
            (cdf2, pdf2),
            (t0, lb, ub),
            NewtonBracket(),
            p=u;
            verbose=false,
            xatol=10e-7,
            atol=10e-7
        )

        lb = ran
        xrand[i] = ran
    end

    shuffle(xrand)
end



module Hansen

    using FFTW
    using LinearAlgebra

    export get_hansen

    # Precision types
    prec = BigFloat
    precc = Complex{BigFloat}

    """
        get_hansen(ecc)

    Wrapper used in your original code:
    computes Hansen coefficients for given eccentricity.
    Relies on global n, m, k_min, k_max.
    """
    function get_hansen(ecc, n::Int, m::Int, k_min::Int, k_max::Int)::Tuple{Array{Float64,1}, Array{Float64,1}}
        k_range2, X = hansen_fft(-(n+1), m, ecc, k_min, k_max; N=2^18)
        return k_range2, X
    end

    """
        nextpow2_int(x)

    Return the integer p such that 2^p ≥ x.
    """
    function nextpow2_int(x::Int)
        return x > 0 ? ceil(Int, log2(x)) : 0
    end


    """
        kepler_newton(M, e)

    Solve Kepler’s equation  E - e*sin(E) = M
    for the eccentric anomaly E, using Newton iteration.
    M may be scalar or array.
    """
    function kepler_newton(M, e)
        M = Float64.(M)
        E = copy(M)

        # Danby-style improved initial guess
        if e > 0
            E .= M .+ e .* sin.(M) ./ (1 .- sin.(M .+ e) .+ sin.(M))
        end

        # Newton iterations
        for _ in 1:10
            f  = E .- e .* sin.(E) .- M
            fp = 1 .- e .* cos.(E)
            dE = -f ./ fp
            E .+= dE
            if maximum(abs.(dE)) < 1e-13
                break
            end
        end

        return mod.(E, 2π)
    end

    """
        hansen_fft(n, m, e, kmin, kmax; N=nothing)

    Compute Hansen coefficients X_k^{n,m}(e) using FFT on mean anomaly.

    Returns a tuple (k, Xkm), where:
    - k is a vector of integer k-indices
    - Xkm are the corresponding Hansen coefficients (real part)
    """
    function hansen_fft(n, m, e, kmin, kmax; N=nothing)

        # Choose FFT size adaptively
        if N === nothing
            width  = max(64, 4*(kmax - kmin + 1))
            target = width * max(8, ceil(Int, 16 / (1 - e + eps())))
            p = max(12, ceil(Int, log2(target)))
            N = 2^p
        else
            p = nextpow2_int(N)
            N = 2^p
        end

        # Mean anomaly grid
        M = (0:N-1) .* (2π / N)

        # Solve Kepler
        E = kepler_newton(M, e)

        ce = cos.(E)
        se = sin.(E)
        r_over_a = 1 .- e .* ce
        v = atan.(sqrt(1 - e^2) .* se, ce .- e)   # true anomaly

        # Hansen integrand
        f = (r_over_a .^ n) .* exp.(im * m .* v)

        # FFT, normalized like Python’s `fft(f)/N`
        F = fftshift(fft(ComplexF64.(f))) ./ N

        k_all = collect(-N÷2 : N÷2-1)
        mask = (k_all .>= kmin) .& (k_all .<= kmax)

        k   = k_all[mask]
        Zk  = F[mask]
        Xkm = real.(Zk)

        return k, Xkm
    end

end
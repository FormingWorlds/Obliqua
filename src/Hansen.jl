# Hansen coefficients X_s^{n,m}(e) of the tidal forcing on an eccentric orbit, by FFT over the mean anomaly.
#
# Kepler's equation is solved with a bracketed Newton iteration (converges for every e < 1), the FFT length follows
# from the eccentricity and is checked by the spectral energy left in the upper half of the band (aliasing), and the
# forcing modes are selected by the share of the spectral weight they hold, so no table of mode ranges is needed.

module Hansen

    using FFTW

    export get_hansen, hansen_modes, hansen_spectrum, kepler_safe


    """
        nextpow2_int(x)

    Return the integer p such that 2^p ≥ x.
    """
    function nextpow2_int(x::Int)
        return x > 0 ? ceil(Int, log2(x)) : 0
    end


    """
        kepler_safe(M, e)

    Eccentric anomaly E in [0, 2π) with E - e sin E = M for M in [0, 2π), 0 <= e < 1. Newton iteration kept inside
    the bracket [M, M + e] (for M <= π; the other half by symmetry), with bisection whenever a step leaves the
    bracket, so the iteration converges for every e < 1.
    """
    function kepler_safe(M::Float64, e::Float64)::Float64
        m = mod(M, 2π)
        # E(2π - M) = 2π - E(M): solve on [0, π] only
        flip = m > π
        flip && (m = 2π - m)
        # E - e sin E is increasing and E - M lies in [0, e] on [0, π]
        lo, hi = m, min(m + e, Float64(π))
        E = e < 0.8 ? m + e * sin(m) : hi
        E = clamp(E, lo, hi)
        for _ in 1:200
            f = E - e * sin(E) - m
            f > 0 ? (hi = E) : (lo = E)                      # shrink the bracket around the root
            fp = 1 - e * cos(E)
            En = E - f / fp
            (En > lo && En < hi) || (En = 0.5 * (lo + hi))   # bisect when Newton leaves the bracket
            abs(En - E) <= 4eps(Float64) * max(1.0, abs(En)) && (E = En; break)
            E = En
        end
        return flip ? 2π - E : E
    end

    kepler_safe(M::Real, e::Real) = kepler_safe(Float64(M), Float64(e))


    """
        hansen_spectrum(n, m, e; N=nothing, alias_tol=1e-13)

    Hansen coefficients X_s^{n,m}(e) = (1/2π) ∫ (r/a)^n exp(i m v) exp(-i s M) dM for all s in -N/2 : N/2-1, by FFT
    over the mean anomaly M with `kepler_safe`. Without `N` the length starts at 2^p >= 64 (1 - e)^(-3/2) (at least
    2^12) and doubles until the share of the spectral energy in |s| > N/4 is below `alias_tol` (at most 2^26 points).
    Returns (s, X, N).
    """
    function hansen_spectrum(n::Int, m::Int, e::Float64; N::Union{Nothing,Int} = nothing, alias_tol::Float64 = 1e-13)
        0 <= e < 1 || throw(ArgumentError("hansen_spectrum: need 0 <= e < 1, got $e"))
        # the spectrum widens like (1 - e)^(-3/2) (pericentre passage)
        p = N === nothing ? clamp(ceil(Int, log2(64 / (1 - e)^1.5)), 12, 26) : nextpow2_int(N)
        while true
            Np = 2^p
            f = Vector{ComplexF64}(undef, Np)
            q = sqrt(1 - e^2)
            @inbounds for j in 1:Np
                E = kepler_safe((j - 1) * (2π / Np), e)
                sE, cE = sincos(E)
                ra = 1 - e * cE                              # r / a
                # exp(i v) = ((cos E - e) + i sqrt(1 - e^2) sin E) / (r/a)
                ph = complex(cE - e, q * sE) / ra
                f[j] = ra^n * ph^m
            end
            F = fftshift(fft(f)) ./ Np
            X = real.(F)                                     # real for the symmetric integrand
            s = collect(-Np÷2:Np÷2-1)
            # aliasing check: energy in the upper half of the band
            tot = sum(abs2, X)
            hi = sum(abs2(X[j]) for j in eachindex(X) if abs(s[j]) > Np ÷ 4)
            if N !== nothing || hi <= alias_tol * tot || p >= 26
                hi <= alias_tol * tot || @warn "hansen_spectrum: spectrum not resolved with N = 2^$p at e = $e (aliased share $(hi / tot))"
                return s, X, Np
            end
            p += 1
        end
    end


    """
        get_hansen(e, n, m, s_min, s_max)

    Hansen coefficients X_s^{-(n+1),m}(e) of the degree-n, order-m tide for s = s_min:s_max (from `hansen_spectrum`,
    with the FFT length raised if the range does not fit in its band). Returns (s, X).
    """
    function get_hansen(e::Real, n::Int, m::Int, s_min::Int, s_max::Int)::Tuple{Vector{Int}, Vector{Float64}}
        s_max >= s_min || throw(ArgumentError("get_hansen: need s_max >= s_min, got $s_min, $s_max"))
        s, X, Np = hansen_spectrum(-(n + 1), m, Float64(e))
        smax = max(abs(s_min), abs(s_max))
        if smax >= Np ÷ 4
            # range reaches the band edge (only for ranges far wider than the spectrum): longer FFT
            s, X, Np = hansen_spectrum(-(n + 1), m, Float64(e); N = 4 * smax)
        end
        i0 = findfirst(==(s_min), s)
        return collect(s_min:s_max), X[i0:i0+(s_max-s_min)]
    end


    """
        hansen_modes(e, n, m; tail_tol=1e-6)

    Mode range of the degree-n, order-m tide at eccentricity e, with coefficients X_s^{-(n+1),m}(e) from
    `hansen_spectrum`. The contiguous range [s_min, s_max] is the smallest one that leaves less than `tail_tol` of the
    weight sum_s (1 + |s|) X_s^2 outside on each side: X_s^2 weights the tidal potential and |s| X_s^2 the dissipation
    of a response with a frequency-independent phase lag. For m = 0 only s >= 0 is returned (X_-s = X_s; the caller
    counts the mirrored modes twice).

    Returns (s, X) with s = s_min:s_max.
    """
    function hansen_modes(e::Real, n::Int, m::Int; tail_tol::Real = 1e-6)
        0 < tail_tol < 1 || throw(ArgumentError("hansen_modes: need 0 < tail_tol < 1, got $tail_tol"))
        s, X, _ = hansen_spectrum(-(n + 1), m, Float64(e))
        w = (1 .+ abs.(s)) .* abs2.(X)
        tot = sum(w)
        # drop modes from each end while the dropped weight stays below tail_tol of the total
        acc = 0.0; i0 = 1
        while i0 < length(w) && acc + w[i0] < tail_tol * tot
            acc += w[i0]; i0 += 1
        end
        acc = 0.0; i1 = length(w)
        while i1 > i0 && acc + w[i1] < tail_tol * tot
            acc += w[i1]; i1 -= 1
        end
        if m == 0
            # mirrored spectrum: keep s >= 0 up to the larger of the two ends
            smax = max(abs(s[i0]), abs(s[i1]))
            i0 = findfirst(==(0), s); i1 = findfirst(==(smax), s)
        end
        return s[i0:i1], X[i0:i1]
    end

end

using Test
using Obliqua

# Tests of the eccentricity spectrum (src/Hansen.jl): the bracketed Kepler solver, the FFT spectrum with automatic
# length (hansen_spectrum) against closed forms and an independent quadrature, the coefficients of a given mode range
# (get_hansen) and the mode selection by spectral weight (hansen_modes).

const HN = Obliqua.Hansen

# X_s^{n,m}(e) by the trapezoidal rule over the eccentric anomaly E (dM = (1 - e cos E) dE, M = E - e sin E): no
# Kepler solve, spectrally accurate for a periodic integrand
function hansen_quad(n, m, e, s; N = 2^16)
    acc = 0.0
    for j in 0:N-1
        E = 2π * j / N
        ra = 1 - e * cos(E)
        v = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
        acc += ra^n * cos(m * v - s * (E - e * sin(E))) * ra
    end
    return acc / N
end

@testset "Hansen" begin

    @testset "Kepler's equation for any e < 1" begin
        for e in (0.0, 0.3, 0.9, 0.99, 0.9999, 1 - 1e-9), M in (0.0, 1e-9, 1e-3, 0.5, 3.0, π, 4.0, 2π - 1e-6)
            E = HN.kepler_safe(M, e)
            @test 0 <= E <= 2π
            # residual of the equation itself, wrapped to (-π, π]
            res = mod(E - e * sin(E) - M + π, 2π) - π
            @test abs(res) < 5e-15 * max(1.0, M)
        end
        # circular orbit: E = M exactly
        @test HN.kepler_safe(1.234, 0.0) == 1.234
    end

    @testset "spectrum against closed forms" begin
        # sum over s of (X_s^{-3,m})^2 is the orbit average of (a/r)^6 = (1 + 3e^2 + 3e^4/8) / (1 - e^2)^(9/2)
        for e in (0.05, 0.5, 0.9, 0.95, 0.99), m in (0, 2)
            s, X, N = HN.hansen_spectrum(-3, m, e)
            exact = (1 + 3e^2 + 3e^4 / 8) / (1 - e^2)^4.5
            @test isapprox(sum(abs2, X), exact; rtol = 1e-10)
            # X_0^{-3,0} = (1 - e^2)^(-3/2) and X_0^{-3,2} = 0, exact for every e
            X0 = X[findfirst(==(0), s)]
            m == 0 ? (@test isapprox(X0, (1 - e^2)^-1.5; rtol = 1e-10)) : (@test abs(X0) < 1e-10 * sqrt(exact))
            # m = 0: the spectrum is even in s
            m == 0 && @test isapprox(X[findfirst(==(3), s)], X[findfirst(==(-3), s)]; rtol = 1e-9, atol = 1e-12)
        end
        # small-e series of the three leading semidiurnal coefficients (Kaula's G functions)
        e = 0.02
        s, X, _ = HN.hansen_spectrum(-3, 2, e)
        at(j) = X[findfirst(==(j), s)]
        @test isapprox(at(2), 1 - 5e^2 / 2 + 13e^4 / 16; atol = 1e-9)
        @test isapprox(at(1), -e / 2 + e^3 / 16; atol = 1e-9)
        @test isapprox(at(3), 7e / 2 - 123e^3 / 16; atol = 1e-7)      # next term of the series is O(e^5) ~ 1e-8
        # the series are distinguishable: swapping the s = 1 and s = 3 terms is far outside the tolerance
        @test abs(at(1) - (7e / 2 - 123e^3 / 16)) > 1e-2
        # circular orbit: a single line
        s0, X0, _ = HN.hansen_spectrum(-3, 2, 0.0)
        @test X0[findfirst(==(2), s0)] ≈ 1
        @test sum(abs2, X0) ≈ 1
        @test_throws ArgumentError HN.hansen_spectrum(-3, 2, 1.0)
        @test_throws ArgumentError HN.hansen_spectrum(-3, 2, -0.1)
    end

    @testset "spectrum against quadrature over the eccentric anomaly" begin
        # independent route (no Kepler solve), at a moderate and a high eccentricity, across the peak and the tails
        for (e, ss) in ((0.6, (-5, 1, 2, 7, 30)), (0.95, (-20, 2, 40, 200)))
            s, X, _ = HN.hansen_spectrum(-3, 2, e)
            for sj in ss
                Xq = hansen_quad(-3, 2, e, sj; N = 2^18)
                @test isapprox(X[findfirst(==(sj), s)], Xq; atol = 1e-12 * maximum(abs, X))
            end
        end
    end

    @testset "coefficients of a given mode range" begin
        s6, X6, _ = HN.hansen_spectrum(-3, 2, 0.6)
        s, X = HN.get_hansen(0.6, 2, 2, -30, 60)
        @test s == collect(-30:60) && length(X) == 91
        i0 = findfirst(==(-30), s6)
        @test X == X6[i0:i0+90]
        # a range reaching past the resolved band takes a longer FFT; the far tail is zero to rounding
        sw, Xw = HN.get_hansen(0.1, 2, 2, -5000, 5000)
        @test length(Xw) == 10001
        @test maximum(abs, Xw[abs.(sw) .> 1000]) < 1e-14
        @test isapprox(sum(abs2, Xw), (1 + 3 * 0.1^2 + 3 * 0.1^4 / 8) / (1 - 0.1^2)^4.5; rtol = 1e-10)
        # a single mode, as the legacy and full spectra use
        s1, X1 = HN.get_hansen(0.02, 2, 2, 2, 2)
        @test s1 == [2] && isapprox(X1[1], 1 - 5 * 0.02^2 / 2 + 13 * 0.02^4 / 16; atol = 1e-9)
        @test_throws ArgumentError HN.get_hansen(0.3, 2, 2, 3, 1)
    end

    @testset "mode selection by spectral weight" begin
        for e in (0.1, 0.6, 0.95), m in (0, 2), tol in (1e-4, 1e-8)
            s, X = HN.hansen_modes(e, 2, m; tail_tol = tol)
            sa, Xa, _ = HN.hansen_spectrum(-3, m, e)
            w = (1 .+ abs.(sa)) .* abs2.(Xa)
            inside = m == 0 ? (abs.(sa) .<= s[end]) : (sa .>= s[1]) .& (sa .<= s[end])
            @test sum(w[.!inside]) <= 2 * tol * sum(w)      # at most tol on each side
            @test all(diff(s) .== 1)
            m == 0 && @test s[1] == 0
            @test length(X) == length(s)
            # the coefficients are those of the spectrum
            @test X == Xa[findfirst(==(s[1]), sa):findfirst(==(s[end]), sa)]
        end
        # tighter tolerance and larger eccentricity both widen the range
        n(e, tol) = length(HN.hansen_modes(e, 2, 2; tail_tol = tol)[1])
        @test n(0.6, 1e-8) > n(0.6, 1e-4)
        @test n(0.95, 1e-6) > 5 * n(0.6, 1e-6)
        # circular orbit: only s = 2 (m = 2) or s = 0 (m = 0)
        @test HN.hansen_modes(0.0, 2, 2)[1] == [2]
        @test HN.hansen_modes(0.0, 2, 0)[1] == [0]
        @test_throws ArgumentError HN.hansen_modes(0.3, 2, 2; tail_tol = 0.0)
        @test_throws ArgumentError HN.hansen_modes(0.3, 2, 2; tail_tol = 1.0)
    end
end

using Test
using Obliqua
using SpecialFunctions: gamma

# Tests of the rheologies in Obliqua.complex_modulus: closed forms of every model, the limits that connect them
# (Burgers -> Maxwell and Sundberg-Cooper -> Andrade for Δ -> 0, the relaxed plateau x / (1 + Δ) of Burgers, the
# unrelaxed and relaxed ends), the Debye peak of the Kelvin-Voigt element in the compliance, passivity
# (Im x̃ >= 0 for σ > 0), the reality condition x̃(-σ) = conj(x̃(σ)), and the argument checks.

const CM = Obliqua.complex_modulus
const PC = Obliqua.precc

# one layer with modulus x [Pa] and viscosity η [Pa s]: τ_M = η / x = 1e4 s, so σ τ_M spans 1e-4 .. 1e6 below
x0, η0 = 5e10, 5e14
cm(σ, mat; kw...) = CM(collect(Float64, σ), PC[x0], Float64[η0]; material=mat, kw...)[1, :]

@testset "rheology" begin

    σs = 10.0 .^ range(-8, 2; length=41)                   # rad/s
    τM = η0 / x0

    @testset "closed forms" begin
        s = 1im .* σs .* τM
        # Maxwell and Andrade (unchanged forms)
        @test cm(σs, "maxwell") ≈ x0 .* s ./ (1 .+ s) rtol = 1e-13
        @test cm(σs, "andrade"; α=0.3) ≈ x0 .* s ./ (1 .+ s .+ gamma(1.3) .* s .^ 0.7) rtol = 1e-13
        # Burgers from its mechanical analogue, independent of the τ_K algebra in the source: spring x, dashpot η and a
        # Kelvin-Voigt element (spring x / Δ parallel to dashpot ηK = r η) in series, so the compliances add
        Δ, r = 0.5, 1e-3
        Jb = 1 / x0 .+ 1 ./ (1im .* σs .* η0) .+ 1 ./ (x0 / Δ .+ 1im .* σs .* (r * η0))
        @test cm(σs, "burgers"; Δ=Δ, ηK_ratio=r) ≈ 1 ./ Jb rtol = 1e-12
        # Sundberg-Cooper: the same plus the Andrade compliance Γ(1+α) (iσ τ_M)^(-α) / x
        Jsc = Jb .+ gamma(1.3) .* s .^ -0.3 ./ x0
        @test cm(σs, "sundberg-cooper"; α=0.3, Δ=Δ, ηK_ratio=r) ≈ 1 ./ Jsc rtol = 1e-12
        # parallel models
        @test cm(σs, "voigt") ≈ x0 .+ 1im .* σs .* η0 rtol = 1e-15
        @test cm(σs, "newton") ≈ 1im .* σs .* η0 rtol = 1e-15
        @test all(cm(σs, "elastic") .== x0) && all(cm(σs, "none") .== x0)
    end

    @testset "limits between models" begin
        # Δ = 0 removes the Kelvin-Voigt element exactly
        @test cm(σs, "burgers"; Δ=0.0) ≈ cm(σs, "maxwell") rtol = 1e-13
        @test cm(σs, "sundberg-cooper"; Δ=0.0, α=0.3) ≈ cm(σs, "andrade"; α=0.3) rtol = 1e-13
        # unrelaxed (σ τ_M = 1e6) and relaxed (σ τ_M = 1e-4) ends of the series models
        for mat in ("maxwell", "andrade", "burgers", "sundberg-cooper")
            hi, lo = cm([1e2], mat)[1], cm([1e-8], mat)[1]
            # Maxwell is off by 1/(σ τ_M) = 1e-6; the Andrade term Γ(1.3) (σ τ_M)^-0.3 still removes ~1.6 %
            @test abs(hi - x0) < (mat == "maxwell" ? 2e-6 : 3e-2) * x0
            @test abs(lo) < 2e-4 * x0
            @test cm([0.0], mat)[1] == 0                   # fully relaxed at zero frequency
        end
        # Voigt is the static elastic solid at σ = 0, Newton has no static strength
        @test cm([0.0], "voigt")[1] == x0 && cm([0.0], "newton")[1] == 0
    end

    @testset "relaxed plateau of Burgers" begin
        # τ_K = ηK Δ / x = 1e-6 * 0.25 * 1e4 s = 2.5e-3 s ≪ 1/σ = 10 s ≪ τ_M = 1e4 s: x̃ -> x / (1 + Δ)
        Δ, r = 0.25, 1e-6
        μp = cm([0.1], "burgers"; Δ=Δ, ηK_ratio=r)[1]
        # corrections are 1/(σ τ_M) = 1e-3 (imaginary) and σ τ_K = 2.5e-4
        @test isapprox(real(μp), x0 / (1 + Δ); rtol = 2e-3)
        # discriminating: the unrelaxed x and the Kelvin-spring x / Δ are both far outside the tolerance
        @test abs(real(μp) - x0) > 0.1 * x0 && abs(real(μp) - x0 / Δ) > x0
        # at the same σ Maxwell is still unrelaxed
        @test isapprox(real(cm([0.1], "maxwell")[1]), x0; rtol = 2e-6)
    end

    @testset "Debye peak of the Kelvin-Voigt element" begin
        # the element adds Δ / (1 + i σ τ_K) to x / x̃: its loss part peaks at σ τ_K = 1 with height Δ / 2
        Δ, r = 0.4, 1e-4
        τK = r * η0 * Δ / x0                               # 0.4 s
        σp = 10.0 .^ range(-2, 2; length=401) ./ τK
        loss = -imag.(x0 ./ cm(σp, "burgers"; Δ=Δ, ηK_ratio=r) .- x0 ./ cm(σp, "maxwell"))
        imax = argmax(loss)
        @test isapprox(σp[imax] * τK, 1.0; rtol = 0.03)    # grid step 1 % in log
        @test isapprox(loss[imax], Δ / 2; rtol = 1e-3)
        # the peak is set by τ_K, not by τ_M (τ_M = 1e4 s would put it 2.5e4 times lower)
        @test σp[imax] * τM > 1e3
    end

    @testset "passivity and reality condition" begin
        σpm = vcat(-reverse(σs), σs)
        for mat in ("maxwell", "andrade", "burgers", "sundberg-cooper", "voigt", "newton"),
            (Δ, r, α) in ((0.2, 0.02, 0.3), (2.0, 1e-5, 0.2))
            μ = cm(σpm, mat; Δ=Δ, ηK_ratio=r, α=α)
            n = length(σs)
            @test all(imag.(μ[n+1:end]) .>= 0)             # dissipation >= 0 for σ > 0
            @test μ[1:n] ≈ conj.(reverse(μ[n+1:end])) rtol = 1e-13
        end
    end

    @testset "edge cases" begin
        # zero unrelaxed modulus (a fluid layer): nothing to relax in the series models
        z = CM([1e-3], PC[0.0], Float64[1e3]; material="burgers")
        @test z[1] == 0
        # infinite viscosity: Burgers and Sundberg-Cooper are elastic, without NaN
        for mat in ("burgers", "sundberg-cooper")
            e = CM([1e-3, -1e-3], PC[x0], Float64[Inf]; material=mat)
            @test all(isfinite, e) && all(e .≈ x0)
        end
        # several layers and frequencies keep the (layer, frequency) layout
        M = CM([1e-3, 1e-2, 1e-1], PC[x0, 2x0], Float64[η0, 10η0]; material="burgers")
        @test size(M) == (2, 3)
        @test M[2, 1] ≈ CM([1e-3], PC[2x0], Float64[10η0]; material="burgers")[1]
    end

    @testset "argument checks" begin
        @test_throws ArgumentError cm(σs, "kelvin")
        @test_throws ArgumentError cm(σs, "andrade"; α=0.0)
        @test_throws ArgumentError cm(σs, "sundberg-cooper"; α=1.0)
        @test_throws ArgumentError cm(σs, "burgers"; Δ=-0.1)
        @test_throws ArgumentError cm(σs, "burgers"; ηK_ratio=0.0)
        # the Kelvin-Voigt parameters are not checked for models that do not use them
        @test cm([1.0], "maxwell"; Δ=-1.0)[1] ≈ cm([1.0], "maxwell")[1]
    end
end

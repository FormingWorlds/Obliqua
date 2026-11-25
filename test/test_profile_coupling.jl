using Test
using LinearAlgebra
using Random
using Obliqua
using Obliqua.constants

# Coupling between the relaxation solver grid and the interior-model grid, and the segment rules that decide which solver
# sees which layer.
#
#   1. cell_average / cell_average_map: the heating density computed on the fine solver grid is handed back on the
#      interior layers as a volume average, so the integrated heating is conserved. (Sampling at layer centres misses
#      dissipation concentrated inside a layer, e.g. a poroelastic boundary layer.)
#   2. nonnegative_heating: negative heating densities (only possible for non-passive moduli or round-off) are set to zero
#      and reported, not flipped in sign.
#   3. Angular quadrature of the heating density: the grid resolution `res` (constants) must integrate degree-2 patterns
#      to ~1 %; 60 degrees does not.
#   4. End-to-end energy balance through run_solid1d_mush_relax: Love-number dissipation
#          P_love = -(2n+1) R omega / (8 pi G) Im(k2)
#      equals the returned layer heating integrated over the layer volumes.
#   5. get_layers: a liquid layer at the surface always stays a fluid segment; the other merge rules are unchanged.
#   6. smooth_complex_modulus! (optional, off by default): contract of the blend when it is enabled.
#   7. biot_bulk_dissipation: volumetric (compaction) dissipation consistent with the constitutive law the two-phase solver
#      integrates (Hay et al. 2025, Eqs. 29c, 33, 34b, 35b): (omega/2) Im(K_d) |tr(eps) + p/K_s|^2. The cross term is absent
#      from their Eq. 39b, which is exact only for incompressible grains (K_s -> infinity).

const PC_RP   = 6.0e6
const PC_RC   = 3.0e6
const PC_RAD  = Float64[PC_RC, 4.0e6, 5.0e6, PC_RP]
const PC_RHO  = Float64[4500.0, 4200.0, 3800.0]
const PC_BULK = ComplexF64[1.3e11, 1.3e11, 1.1e11]
const PC_MCORE = 4/3*pi*PC_RC^3*1.0e4
const PC_OMEGA = 1e-5

"""Three-layer elastic mantle on a liquid core. `porous` = index of a porous layer (phi 0.2, k 1e-5 m^2), or 0;
`lossy` = index of a soft, lossy shear layer (mu = 1e9 (1 + 0.1 i) Pa) between stiff layers, or 0;
`lossy_kd` = give the porous layer a lossy drained bulk modulus (K_d = 2e9 (1 + 0.3 i) Pa, so alpha is complex)."""
function pc_planet(; porous=0, lossy=0, lossy_kd=false)
    phi = zeros(3); μc = ComplexF64[6e10, 6e10, 5e10]; κd = copy(PC_BULK); η = Float64[1e21, 1e21, 1e19]
    if porous > 0
        phi[porous] = 0.2; μc[porous] = 1e9; κd[porous] = lossy_kd ? ComplexF64(2e9, 6e8) : 2e9; η[porous] = 1e15
    end
    if lossy > 0
        μc[lossy] = ComplexF64(1e9, 1e8)
    end
    α = ComplexF64.(phi .> 0) .* (1 .- κd ./ PC_BULK)
    k = Float64.(phi .> 0) .* 1e-5
    return phi, μc, κd, η, α, k
end

"""Run the production driver and return (heating integrated over the input layers, P_love, kT)."""
function pc_driver(; porous=0, lossy=0, lossy_kd=false, dr_min=1000, dr_max=2000)
    phi, μc, κd, η, α, k = pc_planet(porous=porous, lossy=lossy, lossy_kd=lossy_kd)
    prf, _, _, _, kT, _ = Obliqua.run_solid1d_mush_relax(PC_OMEGA, copy(PC_RHO), copy(PC_RAD), zeros(3), η, μc, copy(PC_BULK),
        κd, phi, α, k, PC_RP, PC_MCORE, 1.0e4, ComplexF64(0), ComplexF64(1.3e11);
        dr_min=dr_min, dr_max=dr_max, n=2, m=2, core="liquid", inertial_terms=false, visc_l=1e2, bulk_l=1e9)
    dV = 4π/3 .* diff(PC_RAD.^3)
    P_love = -(2*2 + 1) * PC_RP / (8π*G) * PC_OMEGA * imag(kT)
    return sum(prf .* dV), P_love, kT
end

"""Same planet through the low-level pipeline: fine-grid heating density, its grid, and kT (for the sampling guard)."""
function pc_fine(; porous=0, dr_min=1000, dr_max=2000)
    phi, μc, κd, η, α, k = pc_planet(porous=porous)
    nb = length(PC_RAD)
    r_grid, ρ, η, μc, κs, κl, κd, α, ηl, ϕ, k, g_grid, M_tot = Obliqua.solid1d_mush_relax.resample_profiles(
        copy(PC_RAD), copy(PC_RHO), η, μc, copy(PC_BULK), fill(1e9, nb), κd, α, fill(1e2, nb), phi, k, PC_MCORE, dr_min, dr_max)
    rc = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
    y, _ = Obliqua.solid1d_mush_relax.compute_y(rc, ρ, g_grid, μc, κs, PC_OMEGA, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2,
        1.0e4, ComplexF64(0), ComplexF64(1.3e11), [PC_RP, M_tot, G]; core="liquid", inertial_terms=false, patch=false)
    SG = Obliqua.solid1d_mush_relax.define_spherical_grid(res, 2, 2)
    Eμ, Eκ, El = Obliqua.solid1d_mush_relax.get_heating_profile(y, r_grid, ρ, g_grid, μc, κs, PC_OMEGA, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2, SG)
    return r_grid, Eμ .+ Eκ .+ El, y[5, end] - 1
end

@testset "profile coupling: solver grid <-> interior layers" begin

    @testset "cell_average conserves the volume integral" begin
        # Non-aligned grids, seeded random positive field: the coarse average times the coarse volumes must reproduce the
        # fine volume integral to round-off (the method is exact by construction).
        rng = MersenneTwister(42)
        r_coarse = Float64[3.0e6, 4.0e6, 4.7e6, 5.55e6, 6.0e6]
        r_fine = collect(range(3.0e6, 6.0e6; length=301)) .+ vcat(0.0, 1e3 .* rand(rng, 299), 0.0)   # jittered interior
        f = 1.0 .+ 9.0 .* rand(rng, 300)
        avg = Obliqua.cell_average(r_fine, f, r_coarse)
        Vf = sum(f .* diff(r_fine.^3)); Vc = sum(avg .* diff(r_coarse.^3))
        @test length(avg) == 4
        @test isapprox(Vc, Vf; rtol=1e-12)
        @test all(minimum(f) - 1e-12 .<= avg .<= maximum(f) + 1e-12)       # an average stays within the data range
    end

    @testset "cell_average is exact for a field constant on each coarse layer" begin
        # Each coarse layer split into 7 fine shells carrying the layer value: the averages are the layer values.
        r_coarse = Float64[3.0e6, 4.0e6, 5.0e6, 6.0e6]
        vals = [2.0, 5.0, 11.0]
        r_fine = reduce(vcat, [collect(range(r_coarse[j], r_coarse[j+1]; length=8))[1:end-1] for j in 1:3])
        push!(r_fine, r_coarse[end])
        f = reduce(vcat, [fill(vals[j], 7) for j in 1:3])
        @test Obliqua.cell_average(r_fine, f, r_coarse) ≈ vals rtol=1e-13
        # fine shells extending beyond the coarse range are ignored (uniform field keeps its value)
        r_wide = collect(range(2.9e6, 6.1e6; length=200))
        @test Obliqua.cell_average(r_wide, fill(3.0, 199), r_coarse) ≈ fill(3.0, 3) rtol=1e-12
    end

    @testset "boundary-layer heating: centre sampling misses it, the average does not" begin
        # One layer 4000-5000 km with heating 1 everywhere except a 2 km boundary layer of heating 1e4 at its top.
        # Exact layer average: (1*(V - Vb) + 1e4*Vb)/V. The centre value is 1. Guard: the two differ by > 10x the tolerance.
        r_fine = vcat(collect(range(4.0e6, 4.998e6; length=500)), 5.0e6)
        f = ones(length(r_fine) - 1); f[end] = 1e4
        V = 5.0e6^3 - 4.0e6^3; Vb = 5.0e6^3 - 4.998e6^3
        exact = (V - Vb + 1e4 * Vb) / V
        avg = Obliqua.cell_average(r_fine, f, [4.0e6, 5.0e6])[1]
        @test avg ≈ exact rtol=1e-12
        @test abs(exact - 1.0) > 10 * 1e-12 * exact      # centre sampling (=1) would be wrong
        @test exact > 1.5                                  # and wrong by a large factor here
    end

    @testset "cell_average error contract" begin
        @test_throws ArgumentError Obliqua.cell_average([1.0, 2.0, 3.0], [1.0], [1.0, 3.0])          # length mismatch
        @test_throws ArgumentError Obliqua.cell_average([1.0, 3.0, 2.0], [1.0, 1.0], [1.0, 3.0])     # not increasing
        @test_throws ArgumentError Obliqua.cell_average([1.0, 2.0], [1.0], [1.0])                    # no coarse shell
        m = Obliqua.cell_average_map([1.0, 2.0, 3.0], reshape([1.0, 3.0], 1, 1, 2), [1.0, 3.0])
        @test size(m) == (1, 1, 1)
        @test m[1, 1, 1] ≈ (1.0 * (8 - 1) + 3.0 * (27 - 8)) / (27 - 1) rtol=1e-13
    end

    @testset "nonnegative_heating: clamp and report, never flip" begin
        r = collect(range(1.0, 2.0; length=11))
        E = fill(2.0, 10)
        @test Obliqua.nonnegative_heating(E, r) == E                                   # passive input untouched
        Eround = copy(E); Eround[3] = -1e-12                                           # round-off level: silent clamp
        @test_logs min_level=Base.CoreLogging.Warn Obliqua.nonnegative_heating(Eround, r)
        @test Obliqua.nonnegative_heating(Eround, r)[3] == 0.0
        Ebad = copy(E); Ebad[5] = -1.0                                                 # genuinely negative: warn
        out = @test_logs (:warn, r"Negative tidal heating") Obliqua.nonnegative_heating(Ebad, r)
        @test out[5] == 0.0
        @test out[4] == 2.0
        @test sum(out) < sum(abs.(Ebad))                                               # not flipped to +1
    end

    @testset "angular quadrature at the configured resolution" begin
        # Sphere area and the integral of |Y_22|^2 on the heating grid, against a 5 degree reference.
        quad(g, F) = sum(sin.(g.clats) .* F) * deg2rad(g.res)^2
        G5  = Obliqua.solid1d_mush_relax.define_spherical_grid(5.0, 2, 2)
        Gr  = Obliqua.solid1d_mush_relax.define_spherical_grid(res, 2, 2)
        G60 = Obliqua.solid1d_mush_relax.define_spherical_grid(60.0, 2, 2)
        area(g) = quad(g, ones(length(g.clats), length(g.lons)))
        Y2(g) = quad(g, abs.(g.Y[1, :, :]).^2)
        @test res <= 20.0
        @test isapprox(area(Gr), 4π; rtol=0.015)
        @test isapprox(Y2(Gr), Y2(G5); rtol=0.02)
        # guard: the former 60 degree grid misses these by more than the tolerance
        @test abs(area(G60) / (4π) - 1) > 0.05
        @test abs(Y2(G60) / Y2(G5) - 1) > 0.02
    end

    @testset "energy balance through run_solid1d_mush_relax" begin
        # Porous surface layer: Darcy dissipation sits in boundary layers (delta ~ 10 km) at the top and bottom of a
        # 1000 km layer. The returned layer heating must close the Love-number budget.
        H, P, kT = pc_driver(porous=3)
        @test imag(kT) < 0
        @test H / P ≈ 1.0 atol=0.03
        # Guard: the former centre sampling of the same fine profile captures only a small part of it.
        r_grid, E, kT2 = pc_fine(porous=3)
        @test kT2 ≈ kT rtol=1e-10
        centres = 0.5 .* (PC_RAD[1:end-1] .+ PC_RAD[2:end])
        rc = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
        sampled = [E[argmin(abs.(rc .- c))] for c in centres]
        H_sampled = sum(sampled .* (4π/3 .* diff(PC_RAD.^3)))
        @test abs(H_sampled / P - 1) > 0.3

        # Soft lossy shear layer sandwiched between stiff layers (sharp rigidity jumps on both sides).
        Hs, Ps, kTs = pc_driver(lossy=2)
        @test imag(kTs) < 0
        @test Hs / Ps ≈ 1.0 atol=0.03

        # Porous layer with a lossy drained bulk modulus (complex alpha): compaction dissipation. With the Eq. 39b form of
        # Hay et al. (2025), |tr eps|^2 + |p/K_s|^2 without the cross term, this ratio was 1.094.
        Hk, Pk, kTk = pc_driver(porous=3, lossy_kd=true)
        @test imag(kTk) < 0
        @test Hk / Pk ≈ 1.0 atol=0.01
        @test abs(1.094 - 1.0) > 5 * 0.01                     # guard: the former value is far outside the tolerance
    end

    @testset "biot_bulk_dissipation follows the two-phase constitutive law" begin
        # Purely volumetric strain (no shear) so that only the compaction work remains. Constitutive law of the solver
        # (Hay et al. 2025): mean stress  s = K_d e - alpha p  (Eq. 35b), fluid content  zeta = alpha e + S p  (Eq. 34b with
        # zeta = -phi tr eps_rel), alpha = 1 - K_d/K_s (Eq. 29c), S = phi/K_l + (alpha - phi)/K_s (Eq. 33).
        # Cycle-averaged work rate for exp(i omega t) fields: (omega/2) Im(s e* + p zeta*).
        bbd = Obliqua.solid1d_mush_relax.common.biot_bulk_dissipation
        rng = MersenneTwister(7)
        ω, ϕ, Kl, Ks = 1e-5, 0.2, 1e9, 1.1e11
        Kd = ComplexF64(2e9, 6e8)
        α = 1 - Kd/Ks
        S = ϕ/Kl + (α - ϕ)/Ks
        e = randn(rng, ComplexF64, 4, 5) .* 1e-8
        p = -(1e2 .+ 1e2 .* rand(rng, 4, 5)) .* e .* 1e9      # pore pressure rises under compression (anti-phase with e)
        work = (ω/2) .* imag.((Kd .* e .- α .* p) .* conj.(e) .+ p .* conj.(α .* e .+ S .* p))
        D = bbd(e, p, ω, Kd, α, ϕ, Kl, Ks)
        @test D ≈ work rtol=1e-10
        @test D ≈ (ω/2) .* imag(Kd) .* abs2.(e .+ p ./ Ks) rtol=1e-10
        @test all(D .>= 0)                                       # passive frame (Im K_d > 0) cannot generate energy
        # guard: Eq. 39b (no cross term) over-estimates here by much more than the tolerance
        D39b = (ω/2) .* imag(Kd) .* (abs2.(e) .+ abs2.(p ./ Ks))
        @test sum(D39b) / sum(D) > 1.05
        # limits in which Eq. 39b is exact: non-porous cell (phi = 0, p = 0) and incompressible grains (K_s -> infinity)
        @test bbd(e, zero(p), ω, Kd, α, 0.0, Kl, Ks) ≈ (ω/2) .* imag(Kd) .* abs2.(e) rtol=1e-12
        Ks∞ = 1e30; α∞ = 1 - Kd/Ks∞
        @test bbd(e, p, ω, Kd, α∞, ϕ, Kl, Ks∞) ≈ (ω/2) .* imag(Kd) .* (abs2.(e) .+ abs2.(p ./ Ks∞)) rtol=1e-9
        # elastic frame (real K_d, real alpha): no compaction dissipation at all
        @test maximum(abs.(bbd(e, p, ω, ComplexF64(2e9), 1 - 2e9/Ks, ϕ, Kl, Ks))) < 1e-12 * maximum(D)
    end

    @testset "get_layers: a liquid surface layer stays a fluid segment" begin
        # 12 layers: thin (40 km) layers at the bottom (1), in the interior (6) and at the top (12), 292 km elsewhere.
        # H = 2748 km, min_frac 0.02 -> 55 km threshold, so only layers 1, 6 and 12 are "thin".
        dz = vcat(40e3, fill(292e3, 4), 40e3, fill(292e3, 5), 40e3)
        r = 3.0e6 .+ vcat(0.0, cumsum(dz))
        N = length(r) - 1
        solid, liquid = 1e20, 1e2
        layers(η) = Obliqua.get_layers(r, η, 5e4, 5e4; min_frac=0.02)

        η = fill(solid, N); η[end] = liquid                                 # thin magma ocean on solid: kept
        @test layers(η)[5] == ["solid", "fluid"]
        @test layers(η)[4] == [(1, N - 1), (N, N)]

        # unchanged rules
        η = fill(liquid, N); η[end] = solid                                 # thin solid lid on thick liquid -> merged
        @test layers(η)[5] == ["fluid"]
        η = fill(solid, N); η[1] = liquid                                   # thin melt layer on the core -> merged
        @test layers(η)[5] == ["solid"]
        η = fill(solid, N); η[6] = liquid                                   # thin melt sliver inside solid -> merged
        @test layers(η)[5] == ["solid"]
        η = fill(solid, N); η[end-2:end] .= liquid                          # thick magma ocean -> fluid segment
        @test layers(η)[5] == ["solid", "fluid"]
        @test layers(η)[4] == [(1, N - 3), (N - 2, N)]
    end

    @testset "smooth_complex_modulus! contract (optional blend)" begin
        r = collect(range(3.0e6, 6.0e6; length=9))                          # 8 layers, boundaries
        μ = reshape(ComplexF64[6e10, 6e10, 6e10, 6e10, 6e10, 1e3 + 1e3im, 10 + 0.1im, 10 + 0.1im], 8, 1)
        κ = reshape(fill(ComplexF64(1e11), 8), 8, 1)
        η_smooth = fill(1e20, 8)                                            # no >3-decade drop: untouched
        μ1, κ1 = Obliqua.smooth_complex_modulus!(copy(μ), copy(κ), r, η_smooth)
        @test μ1 == μ
        η_jump = Float64[1e20, 1e20, 1e20, 1e20, 1e20, 1e12, 1e2, 1e2]      # drops of 8 and 10 decades
        μ2, _ = Obliqua.smooth_complex_modulus!(copy(μ), copy(κ), r, η_jump)
        # window = [first jump - 2, last jump + 3] clipped = layers 3..8; its end points keep their own moduli
        @test μ2[3, 1] ≈ μ[3, 1] rtol=1e-12
        @test μ2[8, 1] ≈ μ[8, 1] rtol=1e-12
        @test μ2[1:2, 1] == μ[1:2, 1]
        # interior points are a log-linear blend: |mu| decreases monotonically through the window
        @test all(diff(abs.(μ2[3:8, 1])) .< 0)
        @test μ2[6, 1] != μ[6, 1]
    end
end

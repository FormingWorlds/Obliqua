using Test
using LinearAlgebra
using Obliqua
using Obliqua.solid1d_mush_relax
using Obliqua.solid1d_mush_relax.common
using Obliqua.constants

# Boundary and interface behaviour of the poroviscoelastic relaxation solver (solid1d_mush_relax).
#
# Test planet: elastic three-layer mantle on a liquid core (R = 6000 km, r_c = 3000 km) with ONE porous
# layer (phi = 0.2) placed on the core, in the mid-mantle, or at the surface. All moduli are real, so
# Darcy flow is the only dissipation channel and every result below has a closed energy budget:
#
#     P_love = -(2n+1) R omega / (8 pi G) Im(k2)   must equal   integral of the Darcy heating over volume,
#
# where the heating is evaluated independently from the same y-functions (positive definite).
#
# Resolution matters: the pore pressure forms a diffusion boundary layer of thickness
# delta = sqrt(2 k M / (eta_l omega)) at every impermeable boundary (M = 1/(S + alpha^2/(K_d + 4 mu/3)) ~ 2e9 Pa,
# the effective Biot modulus). With k = 1e-5 m^2, eta_l = 1e2 Pa s and omega = 1e-5 rad/s, delta ~ 6 km, resolved by the 2 km and 1 km
# grids used here. On grids coarser than delta, Im(k2) stays close to its converged value but the Darcy
# heating of the unresolved boundary layer does not, so the energy budget fails; see the "under-resolved" testset.
#
# Boundary conditions for the Darcy flow: R = 0 at the CMB, at the free surface, at the base of a porous
# layer resting on solid, and at the top of a porous layer buried under solid (pore pressure free). At the
# two interfaces the condition holds at the cell face between the last node on one side and the first node
# on the other (the R row of the interface cell), so R at the porous node next to the face is O(dr).

const MI_RP    = 6.0e6
const MI_RC    = 3.0e6
const MI_RAD   = Float64[MI_RC, 4.0e6, 5.0e6, MI_RP]
const MI_RHO   = Float64[4500.0, 4200.0, 3800.0]
const MI_BULK  = ComplexF64[1.3e11, 1.3e11, 1.1e11]
const MI_MCORE = 4/3*pi*MI_RC^3*1.0e4
const MI_OMEGA = 1e-5
const MI_PERM  = 1e-5

"""Solve the tidal problem with one porous layer at `pos` (1 = on the core, 2 = mid, 3 = surface).
Returns (kT, P_love, E_heat, r_centres, y, phi_on_grid)."""
function mi_solve(pos::Int; perm=MI_PERM, ω=MI_OMEGA, dr=2000)
    phi = zeros(3); phi[pos] = 0.2
    μc = ComplexF64[6e10, 6e10, 5e10]; μc[pos] = 1e9
    κd = copy(MI_BULK); κd[pos] = 2e9
    η = Float64[1e21, 1e21, 1e19]; η[pos] = 1e15
    α = ComplexF64.(phi .> 0) .* (1 .- κd ./ MI_BULK)
    k = Float64.(phi .> 0) .* perm
    κl = fill(1e9, 4); ηl = fill(1e2, 4)
    r_grid, ρ, η, μc, κs, κl, κd, α, ηl, ϕ, k, g_grid, M_tot = solid1d_mush_relax.resample_profiles(
        copy(MI_RAD), copy(MI_RHO), η, μc, copy(MI_BULK), κl, κd, α, ηl, phi, k, MI_MCORE, dr, round(Int, 1.01dr))
    rc = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
    y, _ = solid1d_mush_relax.compute_y(rc, ρ, g_grid, μc, κs, ω, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2,
        1.0e4, ComplexF64(0), ComplexF64(1.3e11), [MI_RP, M_tot, G]; core="liquid", inertial_terms=false, patch=false)
    kT = y[5, end] - 1
    SG = solid1d_mush_relax.define_spherical_grid(20.0, 2, 2)
    Eμ, Eκ, El = solid1d_mush_relax.get_heating_profile(y, r_grid, ρ, g_grid, μc, κs, ω, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2, SG)
    dv = 4π/3 .* diff(r_grid.^3)
    P_love = -(2*2 + 1) * MI_RP / (8π*G) * ω * imag(kT)
    return kT, P_love, sum((Eμ .+ Eκ .+ El) .* dv), rc, y, ϕ
end

"""Indices of the first and last porous node on the grid."""
function mi_interfaces(y, ϕ)
    m = findall(ϕ[1:size(y, 2)] .> 0)
    return first(m), last(m)
end

@testset "solid1d_mush_relax: porous boundaries and interfaces" begin

    @testset "CMB constraint for a porous layer on the core" begin
        # B y = 0 must impose zero radial Darcy displacement (R) at the CMB and leave the pore pressure (P)
        # free, for every core type. Ordering of the 8-vector: (U, V, Phi, P, X, Y, Psi, R).
        Y8 = [1,2,5,6,3,7,4,8]
        for core in ("liquid", "solid", "inertial")
            μcore = core == "liquid" ? ComplexF64(0) : ComplexF64(1.0)
            Ic = solid1d_mush_relax.common.get_Ic(1e-4, 0.5, 2.0, 0.6, μcore, ComplexF64(2.0), core, 2; G0=1.0, Y=Y8)
            B  = solid1d_mush_relax.get_core_bc!(1e-4, 0.5, 2.0, 0.6, μcore, ComplexF64(2.0), core, 2; G0=1.0, Y=Y8)
            @test size(B) == (4, 8)
            @test rank(B) == 4
            @test norm(B * Ic) < 1e-10
            N = nullspace(B)                         # admissible CMB states
            @test maximum(abs.(N[Y8[8], :])) < 1e-12 # R = 0 in every admissible state
            @test maximum(abs.(N[Y8[7], :])) > 0.5   # P is not constrained
        end
    end

    @testset "energy balance, surface porous layer (no buried interface)" begin
        kT, PL, E, rc, y, ϕ = mi_solve(3; dr=2000)
        @test imag(kT) < 0                           # passive
        @test PL / E ≈ 1.0 atol=5e-3                 # Love-number dissipation == Darcy heating
        a, b = mi_interfaces(y, ϕ)
        @test iszero(y[8, a-1])                     # no Darcy displacement in the solid
        @test abs(y[8, b]) < 1e-10 * maximum(abs.(y[8, a:b]))   # impermeable free surface (round-off level)
    end

    @testset "energy balance and passivity, buried porous layers (resolved grid)" begin
        for pos in (1, 2)
            k1, PL1, E1, _, _, _ = mi_solve(pos; dr=2000)
            k2, PL2, E2, _, _, _ = mi_solve(pos; dr=1000)
            @test imag(k1) < 0 && imag(k2) < 0
            err1, err2 = abs(PL1/E1 - 1), abs(PL2/E2 - 1)
            @test err2 < 5e-3                        # budget closes on the 1 km grid (observed 1.8e-4, 2.8e-4)
            @test err2 < 0.7 * err1                  # and converges (first order: ~0.5)
            # Im(k2) converged to better than 2 % between the two grids
            @test isapprox(imag(k2), imag(k1); rtol=2e-2)
        end
    end

    @testset "continuity of the elastic state across solid/mush interfaces" begin
        # U, V, total radial stress, tangential stress, potential and potential stress are continuous:
        # the step across an interface must be comparable to the neighbouring grid steps.
        for pos in (1, 2, 3)
            _, _, _, rc, y, ϕ = mi_solve(pos; dr=1000)
            a, b = mi_interfaces(y, ϕ)
            for i in (a - 1, b)                       # lower and upper interface (node i | i+1)
                (i < 2 || i + 2 > size(y, 2)) && continue
                for c in 1:6
                    jump = abs(y[c, i+1] - y[c, i])
                    ref  = max(abs(y[c, i] - y[c, i-1]), abs(y[c, i+2] - y[c, i+1]))
                    @test jump <= 3 * ref + 1e-14 * maximum(abs.(y[c, :]))
                end
            end
        end
    end

    @testset "base of a porous layer on solid: zero Darcy displacement" begin
        _, _, _, rc, y, ϕ = mi_solve(2; dr=1000)
        a, b = mi_interfaces(y, ϕ)
        Rmax = maximum(abs.(y[8, a:b]))
        @test iszero(y[8, a-1])                     # imposed at the last solid node
        @test abs(y[8, a]) < 0.25 * Rmax             # O(dr) at the first porous node
    end

    @testset "top of a buried porous layer (mush -> solid): impermeable" begin
        # No radial Darcy flux into the overlying solid: R = 0 at the cell face above the last porous node b,
        # pore pressure free. R is extrapolated linearly from nodes b-1, b to the face; it vanishes at second
        # order in dr, while R at node b itself is O(dr) (a closure imposing R = 0 at node b fails the first test).
        R_face(y, rc, b) = y[8, b] + (y[8, b] - y[8, b-1]) * (0.5 * (rc[b] + rc[b+1]) - rc[b]) / (rc[b] - rc[b-1])
        for pos in (1, 2)
            _, _, _, rc, y, ϕ = mi_solve(pos; dr=1000)
            _, _, _, rc_f, y_f, ϕ_f = mi_solve(pos; dr=500)
            a, b = mi_interfaces(y, ϕ)
            _, b_f = mi_interfaces(y_f, ϕ_f)
            Rmax = maximum(abs.(y[8, a:b]))
            Pmax = maximum(abs.(y[7, a:b]))
            Rf, Rf_f = abs(R_face(y, rc, b)), abs(R_face(y_f, rc_f, b_f))
            @test 1e-3 * Rmax < abs(y[8, b]) < 0.1 * Rmax   # R at the last porous node is O(dr), not forced to 0
            @test Rf < 5e-3 * Rmax                          # impermeable at the face (observed 9e-4, 2e-3)
            @test Rf_f < 0.5 * Rf                           # and converging (second order: ~0.27)
            @test abs(y[7, b]) > 0.1 * Pmax          # pore pressure not forced to zero (not drained)
            @test iszero(y[8, b+1]) && iszero(y[7, b+1])   # no porous state in the solid
        end
    end

    @testset "buried-layer Darcy dissipation (regression, resolved grid)" begin
        # Pinned Im(k2) of the cell-face impermeable-top closure (dr = 1 km, k = 1e-5 m^2, omega = 1e-5 rad/s).
        # Tolerance 2e-3: the value is deterministic at fixed grid; the grid dependence between 1 km and
        # 0.5 km is 0.1 % (mid-mantle). Discrimination guards: a drained-type top (R + P = 0) gives
        # Im(k2) = -3.336e-3 (on the core) and -7.669e-4 (mid-mantle), 46x and 30x larger; R = 0 imposed at the
        # last porous node (instead of the face) gives -7.269e-5 and -2.541e-5, 0.5 % and 1.1 % off.
        pinned  = Dict(1 => -7.302432e-05, 2 => -2.514217e-05)
        drained = Dict(1 => -3.3359184589e-03, 2 => -7.6689174643e-04)
        at_node = Dict(1 => -7.2686658148e-05, 2 => -2.5414695886e-05)
        for pos in (1, 2)
            kT, PL, E, _, _, _ = mi_solve(pos; dr=1000)
            @test isapprox(imag(kT), pinned[pos]; rtol=2e-3)
            @test abs(drained[pos] / pinned[pos]) > 10       # guard: the closures are distinguishable
            @test abs(at_node[pos] / pinned[pos] - 1) > 2e-3
            @test !isapprox(imag(kT), drained[pos]; rtol=0.5)
        end
        # A surface porous layer has no buried interface and is unaffected by the closure.
        kT, _, _, _, _, _ = mi_solve(3; dr=1000)
        @test isapprox(imag(kT), -1.2897119681e-06; rtol=2e-2)
    end

    @testset "under-resolved poroelastic boundary layer (documents the failure mode)" begin
        # Grid spacing 20 km >> delta ~ 200 m: Im(k2) is already within ~1 % of its converged value, but the
        # Darcy heating of the unresolved boundary layer is not, so the energy budget fails on coarse grids
        # (observed 37 % at 20 km, 8 % at 10 km, 0.3 % at 2 km).
        k_c, PL_c, E_c, _, _, _ = mi_solve(1; perm=1e-8, dr=20000)
        k_f, PL_f, E_f, _, _, _ = mi_solve(1; perm=1e-8, dr=10000)
        @test imag(k_c) < 0 && imag(k_f) < 0              # passive on both grids
        @test abs(imag(k_f) / imag(k_c) - 1) < 2e-2       # Im(k2) insensitive to the unresolved layer
        @test abs(PL_c/E_c - 1) > 0.1                     # heating budget off on the coarse grid
        @test abs(PL_f/E_f - 1) < 0.5 * abs(PL_c/E_c - 1) # and improves with refinement
    end
end


"""Two adjacent porous layers (mid-mantle and surface) with porosities `phi2`, `phi3` over a solid layer on the core.
Returns (P_love / Darcy-and-frame heating, kT, y, phi on the grid, cell faces)."""
function mi_two(phi2, phi3; perm=MI_PERM, ω=MI_OMEGA, dr=2000)
    phi = Float64[0.0, phi2, phi3]
    μc = ComplexF64[6e10, 1e9, 1e9]; κd = ComplexF64[1.3e11, 2e9, 2e9]; η = Float64[1e21, 1e15, 1e15]
    α = ComplexF64.(phi .> 0) .* (1 .- κd ./ MI_BULK)
    k = Float64.(phi .> 0) .* perm
    κl = fill(1e9, 4); ηl = fill(1e2, 4)
    r_grid, ρ, η, μc, κs, κl, κd, α, ηl, ϕ, k, g_grid, M_tot = solid1d_mush_relax.resample_profiles(
        copy(MI_RAD), copy(MI_RHO), η, μc, copy(MI_BULK), κl, κd, α, ηl, phi, k, MI_MCORE, dr, round(Int, 1.01dr))
    rc = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
    y, _ = solid1d_mush_relax.compute_y(rc, ρ, g_grid, μc, κs, ω, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2,
        1.0e4, ComplexF64(0), ComplexF64(1.3e11), [MI_RP, M_tot, G]; core="liquid", inertial_terms=false, patch=false)
    kT = y[5, end] - 1
    SG = solid1d_mush_relax.define_spherical_grid(20.0, 2, 2)
    Eμ, Eκ, El = solid1d_mush_relax.get_heating_profile(y, r_grid, ρ, g_grid, μc, κs, ω, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2, SG)
    P_love = -(2*2 + 1) * MI_RP / (8π*G) * ω * imag(kT)
    return P_love / sum((Eμ .+ Eκ .+ El) .* (4π/3 .* diff(r_grid.^3))), kT, y, ϕ
end

@testset "solid1d_mush_relax: porosity jump between two porous layers" begin
    # With R (relative melt displacement) continuous across the porosity jump, as before the fix, the budget is
    # off by a grid-independent 16 % (0.1 | 0.4) and 5 % (0.4 | 0.1): P_love / E = 0.838 and 0.950 at dr = 2 and 1 km.
    # Melt flux continuity phi- R- = phi+ R+ (flux_continuity!) closes it to the discretisation error.
    for (p2, p3) in ((0.1, 0.4), (0.4, 0.1))
        c2, k2, _, _ = mi_two(p2, p3; dr=2000)
        c1, k1, y, ϕ = mi_two(p2, p3; dr=1000)
        @test abs(c1 - 1) < 2e-3
        @test abs(c2 - 1) < 2e-3
        @test imag(k1) < 0
        # at the porosity face R jumps by the porosity ratio and phi R is continuous. The nodes sit half a cell
        # from the face inside a ~6 km boundary layer, so the node values match to a few per cent (observed 2-3 %);
        # with R continuous the flux would jump by |phi+ - phi-| / phi = 75-300 %
        j = findfirst(i -> ϕ[i] > 0 && ϕ[i+1] > 0 && ϕ[i] != ϕ[i+1], 1:size(y, 2)-1)
        flux_lo, flux_hi = ϕ[j] * y[8, j], ϕ[j+1] * y[8, j+1]
        @test abs(flux_lo - flux_hi) < 5e-2 * abs(flux_hi)
        @test abs(y[8, j] - y[8, j+1]) > 0.5 * abs(y[8, j+1])     # R itself jumps (ratio 4 or 1/4)
    end
    # no porosity change: identical to the former treatment (flux_continuity! is a no-op)
    C = rand(ComplexF64, 4, 8); D = rand(ComplexF64, 4, 8)
    C2, D2 = solid1d_mush_relax.flux_continuity!(copy(C), copy(D), 0.2, 0.2)
    @test C2 == C && D2 == D
    C3, D3 = solid1d_mush_relax.flux_continuity!(copy(C), copy(D), 0.1, 0.4; Y=[1, 2, 5, 6, 3, 7, 4, 8])
    @test C3[4, :] == 0.1 .* C[4, :] && D3[4, :] == 0.4 .* D[4, :]
    @test C3[1:3, :] == C[1:3, :] && D3[1:3, :] == D[1:3, :]
end


@testset "solid1d_mush_relax: non-porous layers heated with their own bulk modulus" begin
    # Solid nodes are solved with the 6x6 matrix and the grain modulus K_s. A lossy drained modulus K_d ≠ K_s on a
    # non-porous layer (the SCS profile gives one below the porosity threshold) must not add bulk heating. All other
    # moduli are real, so with Im(K_d) counted the heating would exceed P_love by the fictitious bulk term.
    function mi_lossy_kd(kdfac)
        phi = Float64[0.0, 0.0, 0.2]
        μc = ComplexF64[6e10, 6e10, 1e9]; κd = copy(MI_BULK); κd[3] = 2e9; κd[1:2] .*= kdfac
        η = Float64[1e21, 1e21, 1e15]
        α = ComplexF64.(phi .> 0) .* (1 .- κd ./ MI_BULK)
        k = Float64.(phi .> 0) .* MI_PERM
        r_grid, ρ, η, μc, κs, κl, κd, α, ηl, ϕ, k, g_grid, M_tot = solid1d_mush_relax.resample_profiles(
            copy(MI_RAD), copy(MI_RHO), η, μc, copy(MI_BULK), fill(1e9, 4), κd, α, fill(1e2, 4), phi, k, MI_MCORE, 2000, 2020)
        rc = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
        y, _ = solid1d_mush_relax.compute_y(rc, ρ, g_grid, μc, κs, MI_OMEGA, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2,
            1.0e4, ComplexF64(0), ComplexF64(1.3e11), [MI_RP, M_tot, G]; core="liquid", inertial_terms=false, patch=false)
        kT = y[5, end] - 1
        SG = solid1d_mush_relax.define_spherical_grid(20.0, 2, 2)
        Eμ, Eκ, El = solid1d_mush_relax.get_heating_profile(y, r_grid, ρ, g_grid, μc, κs, MI_OMEGA, ρ .* ϕ, κl, κd, α, ηl, ϕ, k, 2, SG)
        P_love = -(2*2 + 1) * MI_RP / (8π*G) * MI_OMEGA * imag(kT)
        return P_love / sum((Eμ .+ Eκ .+ El) .* (4π/3 .* diff(r_grid.^3))), kT
    end
    # compute_y itself is given the lossy K_d (it does not apply the porosity threshold); the interface cell then
    # evaluates the solid-side node with the 8x8 matrix and K_d, so the budget fails: this is the input that
    # run_solid1d_mush_relax now corrects (K_d = K_s and alpha = 0 where phi = 0)
    c0, k0 = mi_lossy_kd(1.0)
    c1, k1 = mi_lossy_kd(ComplexF64(0.5, 0.05))
    @test abs(c0 - 1) < 5e-3
    @test abs(c1 - 1) > 0.5            # observed 2.55: K_d of the solid-side interface node enters the solution
    # through the driver the non-porous layers get K_d = K_s, and the lossy K_d has no effect
    function drv(kdfac)
        phi = Float64[0.0, 0.0, 0.2]
        μc = ComplexF64[6e10, 6e10, 1e9]; κd = copy(MI_BULK); κd[3] = 2e9; κd[1:2] .*= kdfac
        α = ComplexF64.(1 .- κd ./ MI_BULK)
        prf, _, _, _, kT, _ = Obliqua.run_solid1d_mush_relax(MI_OMEGA, copy(MI_RHO), copy(MI_RAD), zeros(3),
            Float64[1e21, 1e21, 1e15], μc, copy(MI_BULK), κd, phi, α, Float64[0, 0, MI_PERM], MI_RP, MI_MCORE, 1.0e4,
            ComplexF64(0), ComplexF64(1.3e11); dr_min=2000, dr_max=2020, core="liquid", inertial_terms=false,
            visc_l=1e2, bulk_l=1e9, porosity_thresh=1e-5)
        P_love = -(2*2 + 1) * MI_RP / (8π*G) * MI_OMEGA * imag(kT)
        return P_love / sum(prf .* (4π/3 .* diff(MI_RAD.^3))), kT
    end
    d0, kd0 = drv(1.0)
    d1, kd1 = drv(ComplexF64(0.5, 0.05))
    @test kd1 == kd0
    @test d1 == d0
    @test abs(d1 - 1) < 1e-2
end

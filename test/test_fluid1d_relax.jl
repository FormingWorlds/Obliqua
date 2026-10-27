using Test
using LinearAlgebra
using JSON
using Obliqua
using Obliqua.constants

# Tests of the whole-column model (module_fluid = "fluid1d_relax", src/fluid1d_relax.jl), against
#   1. the interface table of the paper appendix (number and kind of conditions per face);
#   2. analytic solutions: homogeneous viscoelastic sphere (exact I_C(R) solution, Love's incompressible static k2)
#      with second-order convergence in dr, and a uniform fluid layer on a liquid core against the static
#      (equipotential) fluid k2;
#   3. the interface closure: every face condition holds to round-off in the solved states;
#   4. the energy theorem (heating = tidal work) and passivity (Im k <= 0 for w > 0) on the paper stacks;
#   5. the ObliquaOcean prototype (column_solve) on the same inputs (test planet), and passivity and energy closure on
#      PROTEUS snapshots (test/fixtures/fluid1d_relax_reference.json, from test/fixtures/make_fluid1d_relax_reference.jl);
#   6. the heating maps against the analytic angular integrals, the fluid classification (fluid_mask), the
#      sub-layer refinement and the fluid moduli;
#   7. run_tides with module_fluid = "fluid1d_relax", and the bulk viscosity / melt bulk modulus profiles.

const F1R = Obliqua.fluid1d_relax
const CMN1 = Obliqua.fluid1d_relax.common
const Gc1 = Float64(Obliqua.constants.G)
const REF1 = JSON.parsefile(joinpath(@__DIR__, "fixtures", "fluid1d_relax_reference.json"))

maxwell1(μ, η, ω) = μ * (1im * ω * η) / (μ + 1im * ω * η)
reldiff1(a, b) = abs(a - b) / abs(b)

"""Test planet of the paper appendix (R = 6000 km, core 3000 km, uniform dr): layers [(type, km), ...] from the core."""
function column1(layers; ω = 1e-5, dr = 3e3, R = 6000e3, rc = 3000e3, core = "liquid", slip = :noslip, pores = :closed,
                 σR = 1e-4)
    N = round(Int, (R - rc) / dr)
    r_grid = collect(range(rc, R, length = N + 1)); r_c = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
    tops = rc .+ 1e3 .* cumsum(last.(layers))
    types = [first(layers[findfirst(>(x), tops)]) for x in r_c]
    ρ = zeros(N); μ = zeros(ComplexF64, N); κs = zeros(ComplexF64, N); κd = zeros(ComplexF64, N); α = zeros(ComplexF64, N)
    ϕ = zeros(N); k = zeros(N); ηl = fill(1.0, N); κl = fill(1e9, N); σ = zeros(N)
    for i in 1:N
        if types[i] === :solid
            ρ[i] = 4000; μ[i] = maxwell1(5e10, 1e20, ω); κs[i] = κd[i] = 1e11
        elseif types[i] === :mush
            ρ[i] = 3500; μ[i] = maxwell1(1e9, 1e17, ω); κs[i] = 1e11; κd[i] = 5e9; α[i] = 1 - κd[i] / κs[i]
            ϕ[i] = 0.2; k[i] = 1e-8
        else
            ρ[i] = 3000; μ[i] = 10 + 1im * ω * 100; κs[i] = κd[i] = 1e11; σ[i] = σR
        end
    end
    ρ_core = 10000.0
    g, M_tot = CMN1.get_g(r_grid, ρ, 4π / 3 * rc^3 * ρ_core)
    μ_core, κ_core = core in ("liquid", "inertial-liquid") ? (0.1 + 0im, 5e11 + 0im) : (1e11 + 0im, 2e11 + 0im)
    return (σ = Float64(ω), r_grid = r_grid, r_c = r_c, ρ = ρ, g = g, μ = μ, κs = κs, κl = κl, κd = κd, α = α, ηl = ηl,
            ϕ = ϕ, k = k, ρl = ρ .* ϕ, scales = [R, M_tot, Gc1], ρ_core = ρ_core, μ_core = μ_core, κ_core = κ_core,
            core = core, inertial_terms = true, fluid = types .=== :fluid, σR = σ, types = types, slip = slip, pores = pores)
end

"""Uniform column on a core of the same material ("inertial" I_C), or a uniform fluid layer on a liquid core."""
function uniform1(dr; R = 6000e3, rc = 3000e3, ρ = 4000.0, μ = 5e10 + 0im, K = 1e11 + 0im, ω = 1e-5, type = :solid,
                  core = "inertial", ρ_core = ρ, σR = 0.0)
    N = round(Int, (R - rc) / dr)
    r_grid = collect(range(rc, R, length = N + 1)); r_c = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
    ρv = fill(ρ, N)
    g, M_tot = CMN1.get_g(r_grid, ρv, 4π / 3 * rc^3 * ρ_core)
    z = zeros(N)
    μc, κc = core == "liquid" ? (0.1 + 0im, 5e11 + 0im) : (μ, K)
    return (σ = Float64(ω), r_grid = r_grid, r_c = r_c, ρ = ρv, g = g, μ = fill(ComplexF64(μ), N),
            κs = fill(ComplexF64(K), N), κl = fill(1e9, N), κd = fill(ComplexF64(K), N), α = zeros(ComplexF64, N),
            ηl = ones(N), ϕ = z, k = z, ρl = z, scales = [R, M_tot, Gc1], ρ_core = ρ_core, μ_core = μc, κ_core = κc,
            core = core, inertial_terms = true, fluid = fill(type === :fluid, N), σR = fill(σR, N), types = fill(type, N),
            slip = :noslip, pores = :closed)
end

"""Exact k of the homogeneous compressible inertial sphere: y(R) = I_C(R) c with X = Y = 0, Psi = (2n+1)/R."""
function sphere_exact1(cs; n = 2)
    Y8 = [1, 2, 5, 6, 3, 7, 4, 8]; Y6 = [1, 2, 4, 5, 3, 6]
    R0, M0, ω0, ρ0, G0, g0, μ0, S, _ = CMN1.get_scales(cs.scales...; Y = Y8)
    R = cs.r_grid[end]; gR = cs.g[end]
    Ic = CMN1.get_Ic(cs.σ / ω0, R / R0, cs.ρ[end] / ρ0, gR / g0, cs.μ[end] / μ0, cs.κs[end] / μ0, "inertial", n;
                     G0 = G0, Y = Y6)
    Sd = [S[i, i] for i in 1:8]; std_of = [findfirst(==(s), Y8) for s in 1:8]
    Yp = [Sd[std_of[s]] for s in (1, 2, 3, 5, 6, 7)] .* Ic
    c = Yp[4:6, :] \ ComplexF64[0, 0, (2n + 1) / R]
    return sum(Yp[3, :] .* c) - 1
end

"""Static degree-n k of a uniform fluid layer (ρ) on a uniform core (ρc): CMB and surface on the equipotential."""
function fluid_static_k1(rc, R, ρc, ρ; n = 2)
    rj = [rc, R]; drho = [ρc - ρ, ρ]
    M = [4π / 3 * rc^3 * ρc, 4π / 3 * rc^3 * ρc + 4π / 3 * (R^3 - rc^3) * ρ]; g = Gc1 .* M ./ rj .^ 2
    P = [4π * Gc1 * rj[j] * drho[j] / (2n + 1) * (rj[i] <= rj[j] ? (rj[i] / rj[j])^n : (rj[j] / rj[i])^(n + 1))
         for i in 1:2, j in 1:2]
    u = (Diagonal(g) - P) \ ((rj ./ R) .^ n)
    return P[2, :] ⋅ u
end

"""Observed order of convergence between successive grid spacings."""
order1(h, e) = [log(e[i+1] / e[i]) / log(h[i+1] / h[i]) for i in 1:length(h)-1]

stack_layers1(d) = [(Symbol(l[1]), Float64(l[2])) for l in d["layers"]]
cplx1(v) = complex(v[1], v[2])


@testset "fluid1d_relax" begin

    # =========================================================================
    # 1. interface table (tab:interface_conditions)
    # =========================================================================
    @testset "interface conditions" begin
        # number of conditions x = (M- + M+)/2 for every pair of node types
        for lo in (:solid, :mush, :fluid), hi in (:solid, :mush, :fluid)
            lo === hi && continue
            lo === :fluid && hi === :fluid && continue
            up, low = F1R.interface_conditions(lo, hi; slip = :freeslip, pores = :open, ϕlo = 0.2, ϕhi = 0.2)
            Mlo, Mhi = length(F1R.node_slots(lo)), length(F1R.node_slots(hi))
            @test length(up) == Mlo ÷ 2
            @test length(low) == Mhi ÷ 2
            @test length(up) + length(low) == (Mlo + Mhi) ÷ 2
        end
        e(s) = (v = zeros(8); v[s] = 1.0; v)
        z8 = zeros(8)
        # mush -> solid, closed pores: R- = 0 in slot 4 of block row n, everything else continuous
        up, low = F1R.interface_conditions(:mush, :solid)
        @test up[4] == (e(8), z8)
        @test all(up[i] == (e(s), e(s)) for (i, s) in enumerate((1, 2, 3)))
        @test all(low[i] == (e(s), e(s)) for (i, s) in enumerate((5, 6, 7)))
        # solid -> fluid, free slip: Y- = 0 in the V row, Y+ = 0 in the Y row
        up, low = F1R.interface_conditions(:solid, :fluid; slip = :freeslip)
        @test up[2] == (e(6), z8)
        @test low[2] == (z8, e(6))
        @test up[1] == (e(1), e(1))                      # U continuous (no melt)
        # with no slip the same face is plain continuity
        up_ns, low_ns = F1R.interface_conditions(:solid, :fluid; slip = :noslip)
        @test up_ns[2] == (e(2), e(2)) && low_ns[2] == (e(6), e(6))
        # mush -> fluid, open pores: U- + phi R- = U+ and P- + X+ = 0
        up, low = F1R.interface_conditions(:mush, :fluid; slip = :freeslip, pores = :open, ϕlo = 0.25)
        @test up[1] == (e(1) + 0.25 * e(8), e(1))
        @test up[4] == (e(4), -e(5))
        # fluid -> mush, open pores (mirror): U- = U+ + phi R+ and P+ + X- = 0 in slot 8 of block row n+1
        up, low = F1R.interface_conditions(:fluid, :mush; slip = :freeslip, pores = :open, ϕhi = 0.25)
        @test up[1] == (e(1), e(1) + 0.25 * e(8))
        @test low[4] == (-e(5), e(4))
        # solid -> mush is always closed (open pores apply to mush <-> fluid only)
        up, low = F1R.interface_conditions(:solid, :mush; pores = :open)
        @test low[4] == (z8, e(8))
        @test up[1] == (e(1), e(1))
    end

    # =========================================================================
    # 2a. homogeneous sphere: exact solution and second-order convergence
    # =========================================================================
    @testset "homogeneous sphere, second order" begin
        ω = 1e-5
        μ = maxwell1(5e10, 1e20, ω)
        kex = sphere_exact1(uniform1(10e3; μ = μ, ω = ω))
        hs = [1e3, 3e3, 10e3]
        ks = [F1R.column_solve(uniform1(h; μ = μ, ω = ω)).k for h in hs]
        e_re = abs.(real.(ks) .- real(kex)) ./ abs(kex)
        e_im = abs.(imag.(ks) .- imag(kex)) ./ abs(kex)
        @test e_re[1] < 1e-7                              # 3.1e-8 at 1 km in ObliquaOcean (obliqua_grid_order.md)
        @test all(1.8 .< order1(hs, e_re) .< 2.2)
        @test all(1.8 .< order1(hs, e_im) .< 2.2)
        # Obliqua's former convention (boundary rows at the end nodes, gravity of the outer face) is first order:
        # the second-order test above must fail for it
        kl = [F1R.column_solve(uniform1(h; μ = μ, ω = ω); bc_at = :node, g_at = :outer).k for h in hs]
        e_l = abs.(real.(kl) .- real(kex)) ./ abs(kex)
        @test all(order1(hs, e_l) .< 1.3)
        @test e_l[1] > 100 * e_re[1]
    end

    @testset "incompressible static limit: Love's k2" begin
        # exact solution in the limit K/mu -> infinity, w -> 0 against k2 = (3/2) / (1 + 19 mu / (2 rho g R))
        μ = 5e10 + 0im
        cs = uniform1(10e3; μ = μ, K = 5e10 * 1e6 + 0im, ω = 1e-8)
        love = 1.5 / (1 + 19 * real(μ) / (2 * 4000.0 * cs.g[end] * cs.r_grid[end]))
        @test abs(real(sphere_exact1(cs)) / love - 1) < 1e-5
        kn = F1R.column_solve(uniform1(1e3; μ = μ, K = 5e10 * 1e6 + 0im, ω = 1e-8)).k
        @test abs(real(kn) / love - 1) < 1e-5
        # a 10 % error in the rigidity term (the most plausible formula slip) is far outside this tolerance
        @test abs(1.5 / (1 + 1.1 * 19 * real(μ) / (2 * 4000.0 * cs.g[end] * cs.r_grid[end])) / love - 1) > 1e-2
    end

    # =========================================================================
    # 2b. fluid layer on a liquid core: static fluid k2
    # =========================================================================
    @testset "fluid layer, static k2" begin
        ks = fluid_static_k1(3000e3, 6000e3, 10000.0, 3000.0)
        # nearly incompressible fluid at very low frequency, static liquid core as an interface, with drag
        k_of(h) = F1R.column_solve(uniform1(h; ρ = 3000.0, μ = 10 + 1im * 1e-8 * 100, K = 1e18 + 0im, ω = 1e-8,
                                            type = :fluid, core = "liquid", ρ_core = 10000.0, σR = 1e-4)).k
        hs = [1e3, 3e3, 10e3]
        kk = [k_of(h) for h in hs]
        err = abs.(real.(kk) .- ks) ./ ks
        @test err[1] < 2e-5                               # static fluid k2 matched to 2e-5
        @test all(err .< 1e-3)
        @test all(imag.(kk) .<= 0)                        # drag and viscosity only dissipate
        # the homogeneous-planet value 3/2 (the wrong answer if the density jump at the CMB were ignored) is excluded
        @test abs(real(kk[1]) - 1.5) > 0.1
    end

    # =========================================================================
    # 3.-5. paper stacks: residuals, closure, passivity, ObliquaOcean agreement
    # =========================================================================
    @testset "paper stacks" begin
        for d in REF1["stacks"]
            layers = stack_layers1(d)
            cs = column1(layers; dr = d["dr"], slip = Symbol(d["slip"]), pores = Symbol(d["pores"]))
            res = F1R.solve_column(cs)
            sol = res.sol
            k_ref = cplx1(d["k"])
            # agreement with the prototype on the same inputs (same equations: round-off level)
            @test reldiff1(res.k, k_ref) < 1e-10
            @test abs(res.closure - d["closure"]) < 1e-8
            U_ref = cplx1.(d["U"])
            @test maximum(abs.(sol.y[1, :] .- U_ref)) < 1e-9 * maximum(abs.(U_ref))
            # every interface condition holds to round-off
            fr = F1R.face_residuals(sol)
            @test isempty(fr) || maximum(f.residual for f in fr) < 1e-10
            # energy theorem and passivity
            @test abs(res.closure - 1) < (d["dr"] <= 3e3 ? 1e-3 : 2e-3)
            @test res.passive
            # the liquid core as an interface (4x2 solution) fits the solved CMB state
            @test sol.core.misfit < 1e-10
        end
    end

    @testset "core routes agree for a liquid core" begin
        # free slip at a static liquid core: the CMB as an interface (4x2 core) and the left null space of the 6x3 I_C
        # (with its tangential-slip column, and the free pore pressure for a mush) span the same B_1
        for layers in ([(:solid, 3000.0)], [(:mush, 3000.0)], [(:fluid, 3000.0)])
            cs = column1(layers)
            a = F1R.column_solve(cs; core_route = :ic)
            b = F1R.column_solve(cs; core_route = :interface)
            @test reldiff1(a.k, b.k) < 1e-10
        end
    end

    @testset "porosity jump inside a mush: melt flux continuity" begin
        # two mush layers (1500 km each) of porosity 0.2 and ϕ2, k = 1e-6 m^2: with R continuous (the former
        # treatment) the budget is off by +0.4, +1.3, +3.3 % for ϕ2 = 0.25, 0.3, 0.4, independently of dr
        for ϕ2 in (0.25, 0.4), dr in (3e3, 1e3)
            cs = column1([(:mush, 3000.0)]; dr = dr)
            cs = merge(cs, (k = fill(1e-6, length(cs.k)),))
            up = cs.r_c .> 4500e3
            ϕ = copy(cs.ϕ); ϕ[up] .= ϕ2
            cs = merge(cs, (ϕ = ϕ, ρl = cs.ρ .* ϕ))
            res = F1R.solve_column(cs)
            @test abs(res.closure - 1) < 5e-4
            @test res.passive
            fr = F1R.face_residuals(res.sol)
            @test length(fr) == 8                          # one face, x = (8 + 8)/2 conditions
            @test maximum(f.residual for f in fr) < 1e-10
            # face values: phi R continuous, R jumps by the porosity ratio
            f = res.sol.faces[1]
            ym = f.C * res.sol.z[f.i]; yp = -(f.D * res.sol.z[f.i+1])
            @test isapprox(0.2 * ym[8], ϕ2 * yp[8]; rtol = 1e-9)
            @test abs(ym[8] / yp[8] - ϕ2 / 0.2) < 1e-6
            @test isapprox(ym[4], yp[4]; rtol = 1e-9)    # pore pressure continuous
        end
        # boundary-layer diagnostic: delta / dr per mush node, larger on the finer grid
        b3 = F1R.boundary_layer_resolution(column1([(:mush, 3000.0)]; dr = 3e3))
        b1 = F1R.boundary_layer_resolution(column1([(:mush, 3000.0)]; dr = 1e3))
        @test isapprox(b1.min_ratio, 3 * b3.min_ratio; rtol = 1e-6)
        # closed form for the test-planet mush: K_d = 5e9, K_s = 1e11, K_l = 1e9 Pa, phi = 0.2, k = 1e-8 m^2,
        # eta_l = 1 Pa s, omega = 1e-5 rad/s, mu = Maxwell(1e9 Pa, 1e17 Pa s)
        μm = maxwell1(1e9, 1e17, 1e-5); αm = 1 - 5e9 / 1e11
        MB = 1 / (0.2 / 1e9 + (αm - 0.2) / 1e11 + αm^2 / (5e9 + 4μm / 3))
        @test isapprox(b3.min_ratio, sqrt(2 * 1e-8 * abs(MB) / (1.0 * 1e-5)) / 3e3; rtol = 1e-10)
        @test b3.n_under == count(x -> !isnan(x) && x < 2, b3.ratio)
    end

    @testset "non-porous node: heating uses the bulk modulus of its motion matrix" begin
        # a solid node is solved with the 6x6 matrix and K = κs; a lossy κd ≠ κs on it (as the SCS profile gives
        # below the porosity threshold) must not add bulk dissipation. With Im(κd) counted, the budget of this
        # stack exceeds the tidal work by a large factor (most of the dissipation would be the fictitious bulk term).
        cs = column1([(:solid, 2580.0), (:mush, 210.0), (:fluid, 210.0)]; slip = :freeslip, pores = :open)
        sol = cs.types .=== :solid
        κd = copy(cs.κd); κd[sol] .= ComplexF64(5e10, 5e9)
        res0 = F1R.solve_column(cs)
        res = F1R.solve_column(merge(cs, (κd = κd,)))
        @test res.k == res0.k                         # κd of a solid node does not enter the solution
        @test abs(res.closure - 1) < 1e-3
        @test res.closure == res0.closure
    end

    @testset "interface residuals at several dr" begin
        for dr in (1e3, 10e3)
            cs = column1([(:mush, 2400.0), (:fluid, 300.0), (:mush, 300.0)]; dr = dr, slip = :freeslip, pores = :open)
            fr = F1R.face_residuals(F1R.column_solve(cs))
            @test length(fr) == 2 * 7                     # two faces, x = 7 conditions each
            @test maximum(f.residual for f in fr) < 1e-10
        end
    end

    # =========================================================================
    # 5b. PROTEUS snapshots, from the raw interior profiles (column_spec)
    # =========================================================================
    @testset "PROTEUS snapshots: passivity and energy closure" begin
        # the prototype's fluid moduli (mu_reg, K_ad) are no longer used, so its k is not a reference here; the
        # solver itself is compared with the prototype on the paper stacks above
        for d in REF1["snapshots"]
            c = d["cfg"]
            ω = Float64(d["omega"])
            r = Float64.(d["radius"]); ρ = Float64.(d["rho"]); η = Float64.(d["visc"])
            ϕ = Float64.(d["phi"]); k = Float64.(d["perm"])
            shear = ComplexF64.(d["shear"]); bulk = ComplexF64.(d["bulk"]); bulkd = ComplexF64.(d["bulkd"])
            # moduli as run_tides computes them for the two-phase solvers; bulk viscosity eta / (phi + 1e-15)
            ζ = η ./ (ϕ .+ 1e-15)
            μc = Obliqua.complex_modulus([ω], shear, η; material = c["material_mu"], α = Float64(c["alpha"]))[:, 1]
            κc = Obliqua.complex_modulus([ω], bulk, ζ; material = "elastic", α = Float64(c["alpha"]))[:, 1]
            κdc = Obliqua.complex_modulus([ω], bulkd, ζ; material = c["material_k"], α = Float64(c["alpha"]))[:, 1]
            αB = 1.0 .- κdc ./ κc
            sol = c["solid"]
            ρ_core = Float64(c["struct"]["core_density"])
            cs = F1R.column_spec(ω, r, ρ, η, μc, κc, κdc, αB, fill(Float64(sol["bulk_l"]), length(ρ)),
                                 Float64(c["visc_l"]), ϕ, k, 4π / 3 * r[1]^3 * ρ_core, ρ_core,
                                 ComplexF64(c["struct"]["core_shear"]), ComplexF64(c["struct"]["core_bulk"]);
                                 σ_R = Float64(c["sigma_R"]), dr_min = Int(sol["dr_min"]), dr_max = Int(sol["dr_max"]),
                                 core = sol["core"], porosity_thresh = Float64(sol["porosity_thresh"]))
            @test length(cs.r_c) == d["n_nodes"]
            @test cs.core == d["core"]
            # every layer of the prototype's ocean relaxes faster than the tide, so it is fluid here as well
            @test all(cs.fluid_layers[Bool.(d["fluid_layers"])])
            res = F1R.solve_column(cs)
            @test abs(res.closure - 1) < 1e-2
            @test res.passive
        end
    end

    # =========================================================================
    # 6. driver and run_tides
    # =========================================================================
    @testset "run_fluid1d_relax driver" begin
        # three-layer planet: solid, mush, melt on a liquid core, layer boundaries on round radii
        r = Float64[3000e3, 4800e3, 5400e3, 6000e3]
        ρ = Float64[4000, 3500, 3000]; η = Float64[1e20, 1e17, 100.0]
        ϕ = Float64[0.0, 0.2, 1.0]; k = Float64[0.0, 1e-8, 0.0]
        ω = 1e-5
        μ0 = [5e10, 1e9, 1e10]
        mods(w) = (ComplexF64[maxwell1(μ0[i], η[i], w) for i in 1:3], ComplexF64[1e11, 1e11, 1e11],
                   ComplexF64[1e11, 5e9, 1e11])
        μc, κc, κdc = mods(ω); αB = 1.0 .- κdc ./ κc
        κl = fill(1e9, 3)
        R = r[end]
        core = (4π / 3 * r[1]^3 * 1e4, 1e4, 0.1 + 0im, 5e11 + 0im)
        args = (ρ, r, η, μc, κc, κdc, αB, ϕ, k, κl, zeros(3), R, core...)
        kw = (dr_min = 1000, dr_max = 3000, sigma_R = 1e-4, visc_l = 1.0, porosity_thresh = 1e-5)
        prf, Eμ, Eκ, El, Ed, kT, kL = Obliqua.run_fluid1d_relax(ω, args...; kw...)
        # heating on the input layers integrates to the tidal work (closure of the solution itself, no rescaling)
        dv = 4π / 3 .* (r[2:end] .^ 3 .- r[1:end-1] .^ 3)
        P_tide = 5 * R * ω / (8π * Gc1) * (-imag(kT))
        @test abs(sum(prf .* dv) / P_tide - 1) < 2e-3
        @test -imag(kT) > 0 && kL == 0
        @test all(prf .>= 0)
        # the angular integral of the four channel maps is the profile (maps rescaled to the analytic integrals)
        wgt = sin.(Obliqua.fluid1d_relax.define_spherical_grid(Obliqua.constants.res, 2, 2).clats) .*
              deg2rad(Obliqua.constants.res)^2
        ang(M) = dropdims(sum(M .* wgt, dims = (1, 2)), dims = (1, 2))
        @test ang(Eμ) .+ ang(Eκ) .+ ang(El) .+ ang(Ed) ≈ prf rtol = 1e-10
        # on the solver cells: the melt layer is fluid, drag only in fluid cells, Darcy only in mush cells (on the
        # input layers a solver cell that straddles a layer boundary is shared by cell_average)
        cs = F1R.column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, 1.0, ϕ, k, core...; σ_R = 1e-4, dr_min = 1000,
                             dr_max = 3000, porosity_thresh = 1e-5)
        @test cs.fluid_layers == [false, false, true]
        hc = F1R.solve_column(cs).heating
        @test all(hc.drag[cs.types .!== :fluid] .== 0) && all(hc.drag[cs.types .=== :fluid] .> 0)
        @test all(hc.darcy[cs.types .!== :mush] .== 0) && sum(hc.darcy) > 0
        @test Set(cs.types) == Set([:solid, :mush, :fluid])
        # the deepest layer, two cells away from any fluid, has no drag; the fluid layer does
        @test ang(Ed)[1] == 0 && ang(Ed)[3] > 0
        # negative frequency with the conjugate moduli: conjugate Love number and the same heating, solved directly
        μm, κm, κdm = mods(-ω)
        @test μm ≈ conj.(μc)
        outm = Obliqua.run_fluid1d_relax(-ω, ρ, r, η, μm, κm, κdm, 1.0 .- κdm ./ κm, ϕ, k, κl, zeros(3), R, core...; kw...)
        @test outm[6] ≈ conj(kT) rtol = 1e-9
        @test outm[1] ≈ prf rtol = 1e-7
        # zero frequency: real relaxed response, no heating
        out0 = Obliqua.run_fluid1d_relax(0.0, ρ, r, η, mods(1e-14)..., αB, ϕ, k, κl, zeros(3), R, core...; kw...)
        @test imag(out0[6]) == 0
        @test all(out0[1] .== 0)
        # more drag, more drag heating in the fluid layer
        out2 = Obliqua.run_fluid1d_relax(ω, args...; kw..., sigma_R = 2e-4)
        @test ang(out2[5])[3] > ang(Ed)[3]
    end

    @testset "nodal heating: maps against the analytic angular integrals" begin
        # grid quadrature of the strain of strain_on_grid! converges to shear_integral for every component
        n = 2; L = n * (n + 1)
        for (m, amp) in ((2, (0.3 + 0.1im, 0.0, 0.0, 0.0)), (2, (0.0, 0.7im, 0.0, 0.0)), (0, (0.0, 0.0, 1.0, 0.0)),
                         (2, (0.0, 0.0, 0.0, 0.4 - 0.2im)), (1, (0.2, 0.1im, -0.3, 0.5)))
            SG = CMN1.define_spherical_grid(2.0, n, m)
            ϵ = zeros(ComplexF64, size(SG.Y, 2), size(SG.Y, 3), 6)
            CMN1.strain_on_grid!(ϵ, amp..., SG)
            tr = sum(ϵ[:, :, 1:3], dims = 3)[:, :, 1]
            loc = sum(abs2.(ϵ[:, :, 1:3]), dims = 3)[:, :, 1] .+ 2sum(abs2.(ϵ[:, :, 4:6]), dims = 3)[:, :, 1] .-
                  abs2.(tr) ./ 3
            quad = sum(loc .* sin.(SG.clats) .* deg2rad(2.0)^2)
            @test quad ≈ CMN1.shear_integral(amp..., L) rtol = 2e-3
        end
        # the full column: maps of each channel integrate to the analytic profile on a fine grid
        cs = column1([(:solid, 2400.0), (:mush, 300.0), (:fluid, 300.0)]; slip = :freeslip, pores = :open)
        res = F1R.solve_column(cs)
        SG = CMN1.define_spherical_grid(3.0, 2, 2)
        maps = F1R.column_heating_map(cs, res.sol.y, 2, cs.σ, SG; solid = res.heating.solid_method)
        h = res.heating
        for (M, E) in zip(maps, (h.shear, h.bulk, h.darcy, h.drag))
            A = dropdims(sum(M .* sin.(SG.clats) .* deg2rad(3.0)^2, dims = (1, 2)), dims = (1, 2))
            @test sum(abs.(A .- E)) <= 5e-3 * sum(abs.(E)) + 1e-300
        end
        # the 60 degree grid of the former default is far off, so the test above would catch a quadrature at it
        SG60 = CMN1.define_spherical_grid(60.0, 2, 2)
        M60 = F1R.column_heating_map(cs, res.sol.y, 2, cs.σ, SG60; solid = res.heating.solid_method)[1]
        A60 = dropdims(sum(M60 .* sin.(SG60.clats) .* deg2rad(60.0)^2, dims = (1, 2)), dims = (1, 2))
        @test sum(abs.(A60 .- h.shear)) > 5e-2 * sum(abs.(h.shear))
    end

    @testset "fluid_mask: low-viscosity layers connected to the surface, at every frequency" begin
        ηm = [1e9, 1e20, 1e7, 100.0, 100.0]
        @test F1R.fluid_mask(ηm) == [false, false, true, true, true]           # layer 1 is below the cutoff but buried
        @test F1R.fluid_mask(ηm; visc_fluid = 1e6) == [false, false, false, true, true]
        @test F1R.fluid_mask(ηm[1:2]) == [false, false]                       # viscous surface layer: no ocean
        @test F1R.fluid_mask(fill(100.0, 3)) == trues(3)                      # fully molten
        @test F1R.fluid_mask(Float64[]) == Bool[]
        @test F1R.fluid_mask([1e8]) == [false]                                # the cutoff itself is not fluid

        ω = 1e-5
        # planet: solid, viscous mush (omega eta / mu = 1e-8) and melt on a liquid core. A fluid layer with
        # eta ~ 1e10 Pa s closes the energy budget only to 5-8 % (hence the default cutoff visc_fluid = 1e8)
        r = Float64[3000e3, 4800e3, 5400e3, 6000e3]
        ρ = Float64[4000, 3500, 3000]; η = Float64[1e20, 1e7, 100.0]
        ϕ = Float64[0.0, 0.2, 1.0]; k = Float64[0.0, 1e-8, 0.0]
        μc = ComplexF64[maxwell1(5e10, 1e20, ω), maxwell1(1e10, 1e7, ω), maxwell1(1e10, 100.0, ω)]
        κc = ComplexF64[1e11, 1e11, 1e11]; κdc = ComplexF64[1e11, 5e9, 1e11]; αB = 1.0 .- κdc ./ κc
        κl = fill(1e9, 3)
        core = (4π / 3 * r[1]^3 * 1e4, 1e4, 0.1 + 0im, 5e11 + 0im)
        kw = (σ_R = 1e-3, dr_min = 1000, dr_max = 3000, porosity_thresh = 1e-5)
        spec(; o...) = F1R.column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, 1.0, ϕ, k, core...; kw..., o...)
        cs_on = spec()
        cs_off = spec(visc_fluid = 1e6)                      # eta = 1e7 above the cutoff: the mush stays a mush
        @test cs_on.fluid_layers == [false, true, true] && cs_off.fluid_layers == [false, false, true]
        # the fluid mush carries drag; energy theorem and passivity hold, and the Love number changes
        r_on = F1R.solve_column(cs_on); r_off = F1R.solve_column(cs_off)
        @test sum(r_on.heating.drag[cs_on.lay .== 2]) > 0 && sum(r_off.heating.drag[cs_off.lay .== 2]) == 0
        @test abs(r_on.closure - 1) < 2e-3 && r_on.passive
        @test abs(r_on.k - r_off.k) > 1e-3 * abs(r_off.k)
        # the layering is the same at every frequency: at 1e-7 rad/s and at 1e-2 rad/s (moduli at that frequency)
        for w in (1e-7, 1e-2)
            μw = ComplexF64[maxwell1(5e10, 1e20, w), maxwell1(1e10, 1e7, w), maxwell1(1e10, 100.0, w)]
            @test F1R.column_spec(w, r, ρ, η, μw, κc, κdc, αB, κl, 1.0, ϕ, k, core...; kw...).fluid_layers ==
                  cs_on.fluid_layers
        end
        @test_throws ArgumentError spec(visc_fluid = 0.0)
    end

    @testset "refine_layers: interpolated sub-layers and a moving fluid interface" begin
        ω = 1e-4
        r = Float64[3000e3, 4000e3, 5000e3, 6000e3]
        ρ = Float64[4000, 3500, 3000]; η = Float64[1e20, 1e10, 100.0]
        μc = ComplexF64[maxwell1(5e10, 1e20, ω), maxwell1(1e10, 1e10, ω), maxwell1(1e10, 100.0, ω)]
        κc = ComplexF64[2e11, 1e11, 5e10]; κdc = ComplexF64[2e11, 5e9, 5e10]
        κl = Float64[1e11, 5e10, 5e10]; ϕ = Float64[0.0, 0.3, 1.0]; k = Float64[0.0, 1e-8, 0.0]
        ref(n; ηv = η, μv = μc) = F1R.refine_layers(n, r, ρ, ηv, μv, κc, κdc, κl, ϕ, k)

        # one sub-layer is the input
        s1 = ref(1)
        @test s1.r == r && s1.ρ == ρ && s1.η == η && s1.μc == μc

        # odd n_sub: the middle sub-layer of each layer sits at the layer centre and keeps the layer values
        s5 = ref(5)
        @test length(s5.ρ) == 15 && length(s5.r) == 16 && s5.parent == repeat(1:3, inner = 5)
        mid = 3:5:15
        @test s5.ρ[mid] ≈ ρ rtol = 1e-12
        @test s5.η[mid] ≈ η rtol = 1e-12
        @test s5.μc[mid] ≈ μc rtol = 1e-12
        # viscosity is log-linear in radius between two centres (second difference of log eta vanishes)
        lη = log.(s5.η[3:8])
        @test maximum(abs.(diff(diff(lη)))) < 1e-10 * maximum(abs.(lη))
        # flat beyond the first and last centre
        @test all(isapprox.(s5.ρ[1:3], ρ[1]; rtol = 1e-12)) && all(isapprox.(s5.ρ[13:15], ρ[3]; rtol = 1e-12))
        # passive moduli stay passive, a zero imaginary part between two elastic values stays zero, and the
        # sub-layers at -omega are the conjugates of those at +omega
        @test all(imag.(s5.μc) .>= 0) && all(imag.(s5.κc) .== 0)
        @test ref(5; μv = conj.(μc)).μc ≈ conj.(s5.μc) rtol = 1e-12
        @test_throws ArgumentError ref(0)

        # the fluid interface sits where the interpolated viscosity crosses visc_fluid between the two centres
        ηf = 1e5
        fl20 = F1R.fluid_mask(ref(20).η; visc_fluid = ηf)
        rc = 0.5 .* (r[1:end-1] .+ r[2:end])
        r_exp = rc[2] + (log(ηf) - log(η[2])) / (log(η[3]) - log(η[2])) * (rc[3] - rc[2])
        @test abs(ref(20).r[findfirst(fl20)] - r_exp) <= (r[2] - r[1]) / 20
        @test all(fl20[findfirst(fl20):end])                     # one contiguous ocean

        # the ocean thickness follows the mush viscosity continuously: a 0.1 dex change in the mush viscosity moves
        # the interface by a fraction of a layer, never by a whole layer (the jump of the unrefined column)
        function H(ηm)
            ηv = [1e20, ηm, 100.0]; μv = ComplexF64[μc[1], maxwell1(1e10, ηm, ω), μc[3]]
            s = ref(20; ηv = ηv, μv = μv)
            return r[end] - s.r[findfirst(F1R.fluid_mask(s.η; visc_fluid = ηf))]
        end
        Hs = [H(x) for x in 10 .^ range(6, 10, length = 41)]
        @test all(diff(Hs) .<= 0)                              # stiffer mush, shallower interface
        @test Hs[1] - Hs[end] > 0.1 * (r[2] - r[1])            # and it does move
        @test maximum(abs.(diff(Hs))) < 0.25 * (r[2] - r[1])

        # column_spec with sub-layers: the Love number of a smooth column barely changes, the solution is passive
        # and closes, and the fluid nodes follow the sub-layer mask
        core = (4π / 3 * r[1]^3 * 1e4, 1e4, 0.1 + 0im, 5e11 + 0im)
        kw = (σ_R = 1e-3, dr_min = 1000, dr_max = 3000, porosity_thresh = 1e-5, visc_fluid = ηf)
        αB = 1.0 .- κdc ./ κc
        cs1 = F1R.column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, 1.0, ϕ, k, core...; kw...)
        cs4 = F1R.column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, 1.0, ϕ, k, core...; kw..., n_sub = 4)
        @test length(cs4.fluid_layers) == 12
        r1 = F1R.solve_column(cs1); r4 = F1R.solve_column(cs4)
        @test r4.passive && abs(r4.closure - 1) < 2e-3
        @test abs(r4.k - r1.k) < 0.2 * abs(r1.k)
        # interface within one sub-layer (250 km) of the threshold crossing, not at the layer face (5000 km)
        @test abs(cs4.r_grid[findfirst(cs4.fluid)] - r_exp) < 250e3 + 3e3
        @test abs(cs1.r_grid[findfirst(cs1.fluid)] - r[3]) < 3e3
        @test_throws ArgumentError F1R.column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, 1.0, ϕ, k, core...; kw..., n_sub = 0)
    end

    @testset "fluid moduli: the layer's own shear and the undrained bulk modulus" begin
        ω = 1e-4
        r = Float64[3000e3, 4800e3, 5400e3, 6000e3]
        ρ = Float64[4000, 3500, 3000]; η = Float64[1e20, 1e20, 100.0]
        # the fluid layer carries the Maxwell modulus of the melt (G_inf 1e10 Pa) and a melt bulk modulus
        μc = ComplexF64[maxwell1(5e10, 1e20, ω), maxwell1(5e10, 1e20, ω), maxwell1(1e10, 100.0, ω)]
        κc = ComplexF64[2e11, 2e11, 3e11]; αB = zeros(ComplexF64, 3)
        κl = Float64[8e10, 8e10, 8e10]; k = zeros(3)
        core = (4π / 3 * r[1]^3 * 1e4, 1e4, 0.1 + 0im, 5e11 + 0im)
        kw = (σ_R = 1e-3, dr_min = 1000, dr_max = 3000)
        # a melt layer (empty frame, porosity 0.9) gets the Reuss average of melt and grain, a pure melt (porosity 1)
        # the melt modulus; neither is the grain modulus. Shear is the layer's own modulus, porosity and alpha vanish
        κdu = ComplexF64[2e11, 2e11, 0]
        cu = F1R.column_spec(ω, r, ρ, η, μc, κc, κdu, αB, κl, 1.0, [0.0, 0.0, 0.9], k, core...; kw...)
        q = findfirst(cu.fluid)
        @test cu.fluid_layers == [false, false, true]
        @test cu.μ[q] == μc[3]
        @test cu.κs[q] ≈ 1 / (0.9 / κl[3] + 0.1 / real(κc[3])) rtol = 1e-12
        @test cu.κd[q] == cu.κs[q] && cu.α[q] == 0 && cu.ϕ[q] == 0
        @test abs(cu.κs[q] - κc[3]) > 0.5 * real(κc[3])
        c1 = F1R.column_spec(ω, r, ρ, η, μc, κc, κdu, αB, κl, 1.0, [0.0, 0.0, 1.0], k, core...; kw...)
        @test c1.κs[q] ≈ κl[3] rtol = 1e-12
        res_u = F1R.solve_column(cu)
        @test res_u.passive && abs(res_u.closure - 1) < 1e-2
        # the interior's undrained modulus replaces the Gassmann modulus built here, in fluid and non-porous layers,
        # and also in the sub-layers; a non-positive entry falls back to Gassmann
        κu = [1.7e11, 0.0, 6.0e10]
        cv = F1R.column_spec(ω, r, ρ, η, μc, κc, κdu, αB, κl, 1.0, [0.0, 0.0, 0.9], k, core...; kw..., κu = κu,
                             porosity_thresh = 1e-5)
        @test cv.κs[q] == 6.0e10 && cv.κd[q] == 6.0e10 && cv.ϕ[q] == 0
        @test abs(cv.κs[q] - cu.κs[q]) > 0.1 * abs(cu.κs[q])      # not the Gassmann value of cu
        s1 = findfirst(==(1), cv.lay)
        @test cv.κs[s1] == 1.7e11                                  # non-porous layer: the interior's value
        s2 = findfirst(==(2), cv.lay)
        @test cv.κs[s2] == cu.κs[findfirst(==(2), cu.lay)]          # zero entry: own Gassmann (= K_s, no melt)
        cv4 = F1R.column_spec(ω, r, ρ, η, μc, κc, κdu, αB, κl, 1.0, [0.0, 0.0, 0.9], k, core...; kw..., κu = κu, n_sub = 4)
        @test cv4.κs[end] ≈ 6.0e10 rtol = 1e-12                     # flat beyond the last centre
        res_v = F1R.solve_column(cv)
        @test res_v.passive && abs(res_v.closure - 1) < 1e-2
        @test_throws ArgumentError F1R.column_spec(ω, r, ρ, η, μc, κc, κdu, αB, κl, 1.0, [0.0, 0.0, 0.9], k, core...; kw...,
                                                   κu = κu[1:2])
        # a fluid layer needs a shear modulus: zero would make its motion matrix singular
        μ0 = copy(μc); μ0[3] = 0
        @test_throws ArgumentError F1R.column_spec(ω, r, ρ, η, μ0, κc, κdu, αB, κl, 1.0, [0.0, 0.0, 1.0], k, core...; kw...)
    end

    @testset "run_tides: whole column and new profiles" begin
        interior_json_path = joinpath(RES_DIR, "interior_data", "runtests_mantle.json")
        cfg = Obliqua.open_config(joinpath(@__DIR__, "test.toml"))
        cfg["orbit"]["obliqua"]["spectrum"] = "adaptive"
        cfg["orbit"]["obliqua"]["s_min"] = "none"
        cfg["orbit"]["obliqua"]["s_max"] = "none"
        cfg["orbit"]["obliqua"]["module_mushy"] = "none"
        # viscoelastic drained frame, so that the bulk viscosity enters the complex bulk moduli
        cfg["orbit"]["obliqua"]["material_k"] = "andrade"
        omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, phi, ncalc =
            load.load_interior_mush_full(interior_json_path, false)
        perm = Obliqua.interior.get_permeability(phi, cfg)
        perm, phi = Obliqua.interior.limit_porosity(perm, phi, cfg)
        bulkd = Obliqua.interior.get_drained_bulk(bulk, phi, cfg)
        bulk_l = Float64(cfg["orbit"]["obliqua"]["solid"]["bulk_l"])

        # former 14-argument call == the full call with the former profiles (bulk viscosity eta/phi, scalar bulk_l)
        cfg["orbit"]["obliqua"]["module_solid"] = "solid1d-mush-relax"
        cfg["orbit"]["obliqua"]["module_fluid"] = "fluid0d"
        out14 = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, shear, bulk, bulkd, phi, perm, cfg)
        out16 = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, visc ./ (phi .+ 1e-15), shear,
                                  bulk, bulkd, fill(bulk_l, length(rho)), phi, perm, cfg)
        @test out14[5] == out16[5]
        @test out14[1] == out16[1]
        # the melt bulk modulus profile enters the mush solver: a softer melt changes k where a layer is porous
        if any(phi .> Float64(cfg["orbit"]["obliqua"]["solid"]["porosity_thresh"]))
            outKl = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, visc ./ (phi .+ 1e-15), shear,
                                      bulk, bulkd, fill(bulk_l / 10, length(rho)), phi, perm, cfg)
            @test maximum(abs.(outKl[5] .- out16[5])) > 0
        end
        # the bulk viscosity profile enters the complex bulk moduli
        outζ = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, 1e-3 .* visc ./ (phi .+ 1e-15),
                                 shear, bulk, bulkd, fill(bulk_l, length(rho)), phi, perm, cfg)
        @test maximum(abs.(outζ[5] .- out16[5])) > 0
        # a zero bulk viscosity (PROTEUS without an elastic profile) falls back to the former law visc / phi
        out0 = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, zeros(length(rho)), shear, bulk,
                                 bulkd, fill(bulk_l, length(rho)), phi, perm, cfg)
        @test out0[5] == out16[5]
        # the 16-argument call is the 17-argument call with a zero undrained modulus (Obliqua's own Gassmann law)
        out17 = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, visc ./ (phi .+ 1e-15), shear,
                                  bulk, bulkd, fill(bulk_l, length(rho)), zeros(length(rho)), phi, perm, cfg)
        @test out17[5] == out16[5]
        @test_throws ArgumentError Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc,
                                                     visc ./ (phi .+ 1e-15), shear, bulk, bulkd,
                                                     fill(bulk_l, length(rho)), zeros(2), phi, perm, cfg)

        # whole-column model: one fluid segment spanning the interior, layers classified by fluid1d_relax itself
        cfg["orbit"]["obliqua"]["module_fluid"] = "fluid1d_relax"
        # with the segments of test.toml's visc_lus / visc_sus it refuses to run
        @test_throws String Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, visc ./ (phi .+ 1e-15),
                                              shear, bulk, bulkd, fill(bulk_l, length(rho)), phi, perm, cfg)
        cfg["orbit"]["obliqua"]["visc_lus"] = 1e30
        cfg["orbit"]["obliqua"]["visc_sus"] = 1e30
        cfg["orbit"]["obliqua"]["enforce_ec"] = false
        # m = 0 adds modes with negative (and zero) forcing frequency, sigma = -s * omega
        cfg["orbit"]["obliqua"]["m"] = [0, 2]
        # resolve the poroelastic boundary layer at the solid -> mush front near the CMB (with the front an exact grid
        # face, test.toml's 300/3000 m grid leaves a ~4 % closure error there; 100/1000 m: < 0.5 %)
        cfg["orbit"]["obliqua"]["solid"]["dr_min"] = 100
        cfg["orbit"]["obliqua"]["solid"]["dr_max"] = 1000
        run_col() = Obliqua.run_tides(omega, axial, ecc, sma, S_mass, rho, radius, visc, visc ./ (phi .+ 1e-15), shear,
                                      bulk, bulkd, fill(bulk_l, length(rho)), phi, perm, cfg)
        dv = 4π / 3 .* (radius[2:end] .^ 3 .- radius[1:end-1] .^ 3)
        # without rescaling, the heating of the solution itself closes the energy budget to the discretisation error
        prf_raw, blk_raw, _, _, _ = run_col()
        @test abs(sum(prf_raw .* rho .* dv) / blk_raw - 1) < 2e-2
        cfg["orbit"]["obliqua"]["enforce_ec"] = true
        power_prf, power_blk, nmk, σ_range, LNk = run_col()
        @test all(isfinite, LNk)
        @test any(σ_range .< 0) && any(σ_range .> 0)
        # strictly dissipative at every non-zero forcing frequency of either sign (a mode clamped by the final
        # passivity backstop would have Im k = 0 exactly and fail here)
        @test all(sign.(σ_range[.!iszero.(σ_range)]) .* imag.(LNk[.!iszero.(σ_range)]) .< 0)
        @test all(imag.(LNk[iszero.(σ_range)]) .== 0)                         # real k at zero frequency
        @test power_blk > 0
        @test all(power_prf .>= 0)
        # profile and Love-number heating agree (enforce_ec): mass-specific profile x layer mass = total
        @test isapprox(sum(power_prf .* rho .* dv), power_blk; rtol = 1e-8)
    end
end

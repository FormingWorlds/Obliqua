# Non-rotating, radially resolved column model of solid, mush and fluid layers (module_fluid = "fluid1d_relax").
#
# Specification: the paper appendix "Motion matrices and Core solution vector" and "Numerical Strategies"
# (ObliquaOcean docs/appendix_methods.tex). Equation labels in the comments below are the LaTeX labels of that
# appendix, e.g. (eq:A_drag), (eq:interface_LG), (tab:interface_conditions).
#
# One Henyey (Kobayashi 2007) recursion per forcing frequency runs from the core-mantle boundary (CMB) to the free
# surface through every layer of the interior: solid (6 states), mush (8 states) and fluid (6 states), in any order.
# Tilt, loading and self-attraction of a fluid layer are therefore not separate terms; they follow from solving the
# whole column at once. Without spin the spherical-harmonic degrees decouple, so only the forcing degree n is solved,
# and the Love number k is read at the surface R directly (no combination of segments).
#
# Layer physics (paper configuration):
#   solid  : Obliqua's 6x6 viscoelastic motion matrix (eq:motion_matrix6) with the complex moduli of the layer.
#   mush   : Obliqua's 8x8 poro-viscoelastic matrix (eq:motion_matrix8), Hay et al. 2019/2025.
#   fluid  : the 6x6 matrix with the layer's own complex shear modulus (for a melt, mu_c ~ i w eta), the undrained
#            (Gassmann) bulk modulus of its frame, grain and melt moduli (kappa_l for phi = 1), no porosity, and an
#            absolute Rayleigh drag sigma_R on the horizontal motion only (eq:A_drag): rho w^2 -> rho w (w - i sigma_R)
#            in the tangential momentum row. The fluid layers are the surface-connected layers selected by
#            fluid_criterion (fluid_measure): below a viscosity cutoff, or where the drained frame modulus is zero
#            (skeleton failure). With the skeleton criterion the drag also acts on the solid and mush layers below,
#            down to where the loss modulus stops exceeding the storage modulus at the forcing frequency.
#
# Interfaces (sec:general_interface): every face is a list of linear conditions L y- = G y+ on the face values
# y- = C_n y_n and y+ = -D_{n+1} y_{n+1} (eq:face_values, eq:interface_LG), which turns the cell equation into
# (L C_n) y_n + (G D_{n+1}) y_{n+1} = 0 (eq:interface_cell). Paper subset (tab:interface_conditions):
#   free slip at every solid/mush <-> fluid face and at a liquid core (Y = 0 on the side(s) that carry shear, V free),
#   closed pores at solid <-> mush (W = 0 on the porous side, P free),
#   open pores at mush <-> fluid (U_f = U_s + W, P + X_f = 0, with the melt flux W = phi R),
#   melt flux W = phi R continuous where the porosity changes between two mush cells (everything else continuous).
# Core: the CMB rows are the left null space of Obliqua's core solutions I_C (core_rows); a liquid core's I_C carries
# the tangential slip as its own column (free slip). A liquid core is solved as "inertial-liquid", with the static
# "liquid" core where the inertial solution is not finite (|omega| below about 1e-11 rad/s).
#
# Grid (sec "Surface solve and grid", docs/obliqua_grid_order.md): nodes at the cell centres of Obliqua's exponential
# grid; the CMB rows act on the face extrapolation y(r_C) = (I - h_1 A_1) y_1 and the surface rows on
# y(R) = (I + h_N A_N) y_N (eq:lower_bound, eq:upper_bound); gravity is evaluated at the nodes from the enclosed
# mass; k and h are read at R. With these choices k converges at second order in dr.
#
# Heating of the same solution (common.get_heating_profile_nodal): shear, Biot bulk and Darcy, with the
# displacement-based strain in the fluid (y4/mu_f is ill conditioned for mu_f -> 0), and drag 1/2 rho sigma_R w^2 |u_h|^2.
# The energy theorem (2n+1) R w / (8 pi G) (-Im k) = integral of the dissipation is the acceptance check.
#
# State vectors and slots. Physical ("standard") order, as everywhere in Obliqua:
#   y = (y1, ..., y8) = (U, V, sigma_rr, sigma_rtheta, Phi, Psi, P, W)
#   U, V   radial and tangential displacement [m]
#   X = y3 = sigma_rr, Y = y4 = sigma_rtheta [Pa]
#   Phi = y5 potential perturbation, Psi = y6 its gradient term (Obliqua's sign convention, k = y5(R) - 1)
#   P = y7 pore pressure [Pa], W = y8 radial melt flux relative to the matrix, W = phi R with R the relative melt
#   displacement [m] (common.get_A): continuous wherever the pores connect, so a porosity change is no face
# The solver works in the "slot" order of the appendix (eq:y_order), U V Phi P | X Y Psi W, so that the upper half
# (slots 1-4) holds the quantities bounded from below and the lower half (slots 5-8) the tractions:
#   storage slot   1  2  3    4  5  6  7    8
#   component      U  V  Phi  P  X  Y  Psi  W
# A solid or fluid node carries only the six slots [1, 2, 3, 5, 6, 7]; a mush node all eight. The permutation is
# passed to Obliqua's get_A as Y (Y[k] = slot of standard component k): Y6 = [1, 2, 4, 5, 3, 6] for the 6-state
# order (U, V, Phi, X, Y, Psi) and Y8 = [1, 2, 5, 6, 3, 7, 4, 8] for the 8-state order.

module fluid1d_relax

    include("constants.jl")
    include("common.jl")
    include("helper.jl")
    using .constants
    using .common
    using .helper

    using LinearAlgebra
    using Printf

    export column_spec, column_solve, column_heating, column_heating_map, solve_column, face_residuals, tidal_work
    export interface_conditions, interface_step, core_rows, left_null_rows, node_slots, boundary_layer_resolution

    # storage slots: U 1, V 2, Phi 3, P 4, X 5, Y 6, Psi 7, W 8
    const SU, SV, SΦ, SP, SX, SY, SΨ, SR = 1, 2, 3, 4, 5, 6, 7, 8

    # get_A permutations: Y[k] is the slot that standard component k (U, V, y3, y4, y5, y6[, P, W]) occupies
    const Y6 = [1, 2, 4, 5, 3, 6]          # 6-state order U, V, Phi, X, Y, Psi
    const Y8 = [1, 2, 5, 6, 3, 7, 4, 8]    # 8-state order U, V, Phi, P, X, Y, Psi, W


    """
        node_slots(t)

    Storage slots carried by a node of type `t`: all eight for a mush, the six non-porous slots [1, 2, 3, 5, 6, 7]
    (U, V, Phi, X, Y, Psi) for a solid or fluid node. Each node keeps its own size in the recursion.
    """
    node_slots(t::Symbol) = t === :mush ? collect(1:8) : [SU, SV, SΦ, SX, SY, SΨ]


    """Unit vector e_s of C^8 (8-slot storage order)."""
    unit8(s::Int) = (v = zeros(8); v[s] = 1.0; v)


    """
        column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, ηl, ϕ, k, m_core, ρ_core, μ_core, κ_core; kwargs...)

    Frequency-dependent column of the whole interior (core to surface) on Obliqua's exponential grid, as a NamedTuple
    for `column_solve` and `column_heating`. Nodes sit at the cell centres; every node carries its own material.

    The surface-connected nodes selected by `fluid_criterion` are fluid (`fluid_measure`, `fluid_mask`): they keep their
    shear modulus, take the undrained (Gassmann) modulus of their frame, grain and melt moduli as their only bulk modulus
    (the mush limit without Darcy flow; the melt modulus for phi = 1), no porosity, and the drag sigma_R. Porosity below
    `porosity_thresh` is set to zero (as in the mush solver), so a node with phi < porosity_thresh is solid. With `κu`
    (the undrained modulus of the interior rheology, per layer) that modulus replaces the Gassmann modulus built here
    wherever it is positive; the interior computes it with the melt volume fraction.

    # Arguments (per input layer unless noted; r has N_l + 1 boundaries)
    - `ω`                    : Forcing frequency [rad s^-1], nonzero; the moduli are those at ω (conjugate for ω < 0).
    - `r, ρ, η`              : Layer boundaries [m], density [kg m^-3], shear viscosity [Pa s].
    - `μc, κc, κdc, αB`      : Complex shear, grain bulk and drained bulk moduli [Pa] at ω, and the Biot coefficient (the
                               nodes take 1 - κd / κs of their own moduli, which is αB for consistent input).
    - `κl`                   : Melt bulk modulus [Pa] per layer.
    - `ηl`                   : Melt viscosity [Pa s] (scalar), used in the Darcy terms of the mush.
    - `ϕ, k`                 : Porosity [-] and permeability [m^2].
    - `m_core, ρ_core`       : Core mass [kg] and density [kg m^-3].
    - `μ_core, κ_core`       : Core moduli [Pa] (only used by solid/inertial core solutions).

    # Keyword Arguments
    - `σ_R=1e-4`             : Rayleigh drag rate in the fluid layers [s^-1], >= 0.
    - `dr_min, dr_max`       : Grid spacing at the surface and at the CMB [m] (Obliqua's exponential grid).
    - `core="liquid"`        : Core type (`get_Ic`); a liquid core is solved as inertial where that is finite (`core_rows`).
    - `inertial_terms=true`  : Inertia in the motion matrices.
    - `porosity_thresh=0.0`  : Porosity below which a node is non-porous.
    - `slip=:freeslip`       : :freeslip or :noslip at solid/mush <-> fluid faces.
    - `pores=:open`          : :open or :closed at mush <-> fluid faces (solid <-> mush is always closed).
    - `visc_fluid=1e8`       : Viscosity [Pa s] below which the nodes connected to the surface are fluid, at every
                               frequency. Layers with η >~ 1e10 Pa s as fluid break the energy closure of the heating.
    - `fluid_criterion=:viscosity`: `:viscosity` (η < visc_fluid) or `:skeleton` (drained frame modulus zero)
    - `κu=nothing`           : Undrained bulk modulus per input layer [Pa] (entries <= 0: Gassmann of the node here).
    - `profile=:nodes`       : `:nodes` interpolates the layer values between the layer centres to every node
                               (`interpolate_layers`), so the material varies from node to node without steps. `:layers`
                               keeps the values constant per input layer, or per sub-layer with `n_sub` (`refine_layers`);
                               every porosity change is then a step with its own pressure-diffusion layer.
    - `n_sub=1`              : Interpolated sub-layers per input layer for `profile = :layers` (not used with `:nodes`).

    A node below the porosity threshold is solved as non-porous with the undrained (Gassmann) bulk modulus of its
    melt-bearing frame, K_u = K_d + alpha^2 / S (common.gassmann), as both its grain and drained modulus. A poorly
    connected mush with little melt is undrained at tidal frequencies, so the bulk stiffness is continuous where a
    node crosses the threshold; without melt K_u = K_s.
    """
    function column_spec(ω::Real, r::AbstractVector, ρ::AbstractVector, η::AbstractVector, μc::AbstractVector,
                         κc::AbstractVector, κdc::AbstractVector, αB::AbstractVector, κl::AbstractVector, ηl::Real,
                         ϕ::AbstractVector, k::AbstractVector, m_core::Real, ρ_core::Real, μ_core::Number,
                         κ_core::Number; σ_R::Real = 1e-4, dr_min::Integer = 300, dr_max::Integer = 3000,
                         core::String = "liquid", inertial_terms::Bool = true,
                         porosity_thresh::Real = 0.0, slip::Symbol = :freeslip, pores::Symbol = :open,
                         visc_fluid::Real = 1e8, n_sub::Integer = 1, κu = nothing,
                         fluid_criterion::Symbol = :viscosity, profile::Symbol = :nodes)
        fn = "column_spec"
        Nl = length(ρ)
        check_arg(ω != 0, "$fn: ω must be nonzero (a zero-frequency mode is the equilibrium tide)")
        check_arg(visc_fluid > 0, "$fn: visc_fluid = $visc_fluid must be > 0")
        check_arg(σ_R >= 0, "$fn: σ_R = $σ_R must be >= 0 (passivity)")
        check_option(fn, "fluid_criterion", fluid_criterion, (:viscosity, :skeleton))
        check_option(fn, "slip", slip, (:freeslip, :noslip))
        check_option(fn, "pores", pores, (:open, :closed))
        check_option(fn, "profile", profile, (:nodes, :layers))
        check_arg(length(r) == Nl + 1, "$fn: r needs length(ρ) + 1 boundaries")
        check_arg(dr_max > dr_min > 0, "$fn: need dr_max > dr_min > 0, got $dr_min, $dr_max")   # ln(dr_max/dr_min) > 0
        check_arg(n_sub >= 1, "$fn: n_sub = $n_sub must be >= 1")
        check_arg(κu === nothing || length(κu) == Nl, "$fn: κu needs one entry per layer")
        profile === :nodes && n_sub > 1 &&
            @warn "$fn: n_sub = $n_sub is not used with profile = :nodes (material interpolated at every node)" maxlog = 1

        # solver grid (cell faces) and nodes at the cell centres; lay: input layer of each node
        r_grid = solver_grid(Float64(r[1]), Float64(r[end]), Int(dr_min), Int(dr_max))
        r_c = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])
        lay = [clamp(searchsortedfirst(r, x) - 1, 1, Nl) for x in r_c]

        # material at the nodes: interpolated to every node, or constant per (sub-)layer
        if profile === :nodes
            m = interpolate_layers(r_grid, r, ρ, η, μc, κc, κdc, κl, ϕ, k; κu = κu)
        else
            src = n_sub > 1 ? refine_layers(n_sub, r, ρ, η, μc, κc, κdc, κl, ϕ, k; κu = κu) :
                              (r = r, ρ = ρ, η = η, μc = μc, κc = κc, κdc = κdc, κl = κl, ϕ = ϕ, k = k, κu = κu)
            j = [clamp(searchsortedfirst(src.r, x) - 1, 1, length(src.ρ)) for x in r_c]   # (sub-)layer of each node
            m = (ρ = src.ρ[j], η = src.η[j], μc = src.μc[j], κc = src.κc[j], κdc = src.κdc[j], κl = src.κl[j],
                 ϕ = src.ϕ[j], k = src.k[j], κu = src.κu === nothing ? nothing : src.κu[j])
        end
        ρn = Vector{prec}(m.ρ); ηn = Vector{prec}(m.η); μn = Vector{precc}(m.μc); κln = Vector{prec}(m.κl)
        κs = Vector{precc}(m.κc); κd = Vector{precc}(m.κdc); ϕn = Vector{prec}(m.ϕ); kn = Vector{prec}(m.k)
        α = 1 .- κd ./ κs

        # fluid nodes: connected to the surface and selected by fluid_criterion (viscosity cutoff or failed frame);
        # with :skeleton the drag also acts below the fluid, in the surface-connected nodes where loss > storage
        fluid = fluid_mask(fluid_measure(ηn, μn; criterion = fluid_criterion, visc_fluid = visc_fluid,
                                         κdc = κd, κc = κs, ϕ = ϕn))
        dragged = fluid_criterion === :skeleton ? fluid .| fluid_mask(fluid_measure(ηn, μn; criterion = :loss)) : fluid
        check_arg(!any(fluid .& iszero.(μn)), "$fn: a fluid layer has a zero shear modulus (give the melt μ and η)")

        # undrained modulus of node i: the interior's where given, else Gassmann with ϕ as porosity
        undrained(i) = (m.κu !== nothing && m.κu[i] > 0) ? precc(m.κu[i]) : gassmann(κd[i], κs[i], κln[i], ϕn[i])
        # below the porosity threshold (as in run_solid1d_mush_relax): non-porous, the undrained modulus as the only
        # bulk modulus (K_u = K_s without melt), alpha = 0
        nopore = ϕn .< porosity_thresh
        for i in findall(nopore)
            κs[i] = κd[i] = undrained(i); α[i] = 0
        end
        ϕn[nopore] .= 0
        # fluid nodes: the undrained modulus (kappa_l at phi = 1) as grain and drained modulus, no porosity
        for i in findall(fluid)
            κs[i] = κd[i] = undrained(i); α[i] = 0
            ϕn[i] = 0; kn[i] = 0
        end

        # gravity at the outer cell faces from the enclosed mass
        g_grid, M_tot = get_g(r_grid, ρn, prec(m_core))
        # Rayleigh drag on the dragged nodes; node types
        σR = [dragged[q] ? Float64(σ_R) : 0.0 for q in eachindex(r_c)]
        types = [fluid[q] ? :fluid : (ϕn[q] > 0 && kn[q] > 0) ? :mush : :solid for q in eachindex(r_c)]

        return (σ = Float64(ω), r_grid = r_grid, r_c = r_c, ρ = ρn, g = g_grid, μ = μn, κs = κs, κl = κln, κd = κd,
                α = α, ηl = fill(prec(ηl), length(r_c)), ϕ = ϕn, k = kn, ρl = ρn .* ϕn, scales = [r[end], M_tot, G],
                M_tot = M_tot, ρ_core = Float64(ρ_core), μ_core = ComplexF64(μ_core), κ_core = ComplexF64(κ_core),
                core = core, inertial_terms = inertial_terms, fluid = fluid, σR = σR, types = types, lay = lay,
                slip = slip, pores = pores, profile = profile)
    end


    """
        interface_conditions(lo, hi; slip=:freeslip, pores=:open) -> (up, low)

    Conditions l . y- = g . y+ (pairs (l, g) of 8-slot vectors) at the face between a node of type `lo` (below) and
    `hi` (above), with the melt flux W = phi R in slot 8 of a mush node. A condition keeps the row of the component
    whose continuity it replaces:

    - `up`  : the M-/2 rows of block row n, i.e. the upper-half slots of the node below: U, V, Phi and, for a mush
              below, slot 4 (P).
    - `low` : the M+/2 rows of block row n+1, i.e. the lower-half slots of the node above: X, Y, Psi and, for a mush
              above, slot 8 (W).

    In total x = (M- + M+)/2 conditions. A one-sided porous pair (P, W) gives one condition, placed in the row of its
    free member: at the top of a mush (slot 4 of block row n) and at its base (slot 8 of block row n+1). A porosity
    change between two mush nodes is not a face: W, P and all other components are continuous there.
    """
    function interface_conditions(lo::Symbol, hi::Symbol; slip::Symbol = :freeslip, pores::Symbol = :open)
        check_arg(lo !== hi, "interface_conditions: no face between two nodes of type $lo")
        z8 = zeros(8)
        # continuity of slot s: e_s . y- = e_s . y+ (the row of an ordinary cell, L = G = I)
        cont(s) = (unit8(s), unit8(s))
        # initialize the rows of the continuity conditions
        rows = Dict(s => cont(s) for s in (SU, SV, SΦ, SX, SY, SΨ))
        # the 4th and 8th rows of block rows n and n+1 are only present for a mush below or above, respectively
        r4 = r8 = nothing

        # helper to check if a node type is fluid
        isfl(t) = t === :fluid

        # interface: mush | ...
        if lo === :mush                                               # top of a porous layer (slot 4, block row n)
            # ... | fluid (open pores)
            if pores === :open && isfl(hi)
                # open pores: pore pressure = fluid pressure, P- + X+ = 0 (P- = -X+, the fluid's sigma_rr)
                r4 = (unit8(SP), -unit8(SX))
                # the melt leaving the pores enters the fluid: U- + W- = U+ (replaces U continuity)
                rows[SU] = (unit8(SU) + unit8(SR), unit8(SU))
            # ... | solid/fluid (closed pores)
            else
                # closed pores: no Darcy flux through the face, W- = 0 (P- is free)
                r4 = (unit8(SR), z8)
            end
        end
        # interface: ... | mush
        if hi === :mush                                               # base of a porous layer (slot 8, block row n+1)
            # fluid | ... (open pores)
            if pores === :open && isfl(lo)
                # open pores, mirror image: P+ + X- = 0 and U- = U+ + W+
                r8 = (-unit8(SX), unit8(SP))
                rows[SU] = (unit8(SU), unit8(SU) + unit8(SR))
            # solid/fluid | ... (closed pores)
            else
                # closed pores: W+ = 0 (P+ is free)
                r8 = (z8, unit8(SR))
            end
        end
        # interface: solid/mush/fluid | solid/mush/fluid (free slip) (symmetric)
        if slip === :freeslip && isfl(lo) != isfl(hi)
            # free slip at a solid/mush <-> fluid face: zero shear traction on both sides, V may jump.
            # Y- = 0 replaces V continuity (V row of block row n) ...
            rows[SV] = (unit8(SY), z8)
            # ... and Y+ = 0 replaces Y continuity (Y row of block row n+1)
            rows[SY] = (z8, unit8(SY))
        end
        up = [rows[SU], rows[SV], rows[SΦ]]
        length(node_slots(lo)) == 8 && push!(up, r4)                      # 4th row of block row n for a mush below
        low = [rows[SX], rows[SY], rows[SΨ]]
        length(node_slots(hi)) == 8 && push!(low, r8)                     # 4th row of block row n+1 for a mush above
        return up, low
    end


    """
        interface_step(C, D, up, low, slo, shi) -> (Ĉu, D̂u, Ĉl, D̂l)

    Cell equation of a face with conditions (up, low): C^ = L C_n, D^ = G D_{n+1}, with L and G
    restricted to the active slots slo (node below) and shi (node above). The rows of `up` form the part that goes
    into block row n (u) and the rows of `low` the part that goes into block row n+1 (l).
    """
    function interface_step(C, D, up, low, slo, shi)
        # stack the condition vectors (restricted to the active slots of each side) into matrices
        mat(cs, k, sl) = permutedims(reduce(hcat, [c[k][sl] for c in cs]))

        # restrict the condition matrices to the active slots of each side
        Lu, Gu, Ll, Gl = mat(up, 1, slo), mat(up, 2, shi), mat(low, 1, slo), mat(low, 2, shi)
        
        # l . y- = g . y+ with y- = C y_n, y+ = -D y_{n+1}  =>  (l C) y_n + (g D) y_{n+1} = 0
        return Lu * C, Gu * D, Ll * C, Gl * D
    end


    """
        left_null_rows(J) -> Matrix

    Rows spanning the left null space of the full-column-rank matrix J (B J = 0, as `permutedims(nullspace(transpose(J)))`),
    computed from J with every column scaled to unit norm. The scaling leaves the null space unchanged; without it the
    regular core solutions of a small core at high degree, whose norms differ by about (r_C / R)^n, fall below the rank
    tolerance of `nullspace` and the core gets one row too many.
    """
    function left_null_rows(J::AbstractMatrix)
        Js = J ./ transpose([max(norm(c), floatmin(Float64)) for c in eachcol(J)])
        B = permutedims(nullspace(transpose(Js)))
        size(B, 1) == size(J, 1) - size(J, 2) ||
            error("left_null_rows: core solutions of rank $(size(J, 1) - size(B, 1)) instead of $(size(J, 2))")
        return B
    end


    """
        core_rows(core, t1, ω, r, ρ, g, μ, K, n; G0) -> (B, core_used)

    Core rows B_1 at the CMB, B_1 y(r_C) = 0, for a first node of type `t1` (non-dimensional inputs): the left null
    space of the core solutions I_C (`get_Ic`, eq:B1_nullspace) in the node's slots; for a mush on the core the free
    pore pressure column e_P is added (the core columns have W = 0: no melt flux through the CMB, P free). This covers
    the free-slip interface of a liquid core, whose I_C carries the tangential slip as its own column. A liquid core
    is solved as "inertial-liquid", the dynamic solution that tends to the static one as ω -> 0 (agreement to 1e-9 at
    1e-10 to 1e-6 rad/s on the test interior); where that is not finite (|ω| below about 1e-11 rad/s) the static
    "liquid" core is used. Returns the rows and the core type used.
    """
    function core_rows(core::String, t1::Symbol, ω, r, ρ, g, μ, K, n::Int; G0)
        Ic = nothing; used = core
        if core == "liquid"
            Ic = try get_Ic(ω, r, ρ, g, μ, K, "inertial-liquid", n; G0 = G0, Y = Y6) catch; nothing end
            Ic !== nothing && all(isfinite, Ic) ? (used = "inertial-liquid") : (Ic = nothing)
        end
        Ic === nothing && (Ic = get_Ic(ω, r, ρ, g, μ, K, core, n; G0 = G0, Y = Y6))
        if t1 === :mush
            # embed the 6-slot I_C in the 8 slots and add the free pore pressure column
            Ic8 = zeros(ComplexF64, 8, 4); Ic8[node_slots(:solid), 1:3] .= Ic; Ic8[SP, 4] = 1
            Ic = Ic8
        end
        return left_null_rows(Ic), used
    end


    """
        column_solve(cs, n=2; bc_at=:face, g_at=:node)

    Response of the column `cs` (column_spec) to a unit degree-n tidal potential, with every material face closed by
    interface_conditions/interface_step. One Henyey recursion from the CMB to the surface.

    # Keyword Arguments
    - `bc_at=:face`      : CMB rows on y(r_C) = (I - h_1 A_1) y_1, surface rows on y(R) = (I + h_N A_N) y_N, k read at
                           R, material faces at their grid radius (second order). :node is Obliqua's convention
                           (boundary rows at the end nodes, faces at node midpoints; first order).
    - `g_at=:node`       : gravity at each node from the enclosed mass (second order); :outer is Obliqua's g of the
                           outer cell face (first order).

    # Returns
    NamedTuple with `k`, `h` (Love numbers at r_top), `y` (8 x N physical profiles, standard order U, V, y3, y4,
    Phi, Psi, P, W), `z` (non-dimensional node states), `faces` (for face_residuals), `slots`, `types`, `finite`,
    `r_top`, `scale`, `core` (the core type used, `core_rows`), `g_node`, `z_C`, `z_R`.
    """
    function column_solve(cs, n::Int = 2; bc_at::Symbol = :face, g_at::Symbol = :node)
        check_option("column_solve", "bc_at", bc_at, (:face, :node))
        check_option("column_solve", "g_at", g_at, (:node, :outer))

        # non-dimensionalization and node quantities
        t = cs.types
        N = length(t)
        check_arg(N >= 2, "column_solve: the column needs at least two nodes")

        # non-dimensionalization (common.get_scales): R0 = planet radius, M0 = planet mass, G0 = G;
        # S is the diagonal state scale, y_physical = S z (8-state order, used for 6-state nodes through their slots)
        R0, M0, ω0, ρ0, G0, g0, μ0, S, Sinv = get_scales(cs.scales[1], cs.scales[2], cs.scales[3]; Y = Y8)
        
        # check grid
        rg = cs.r_grid
        check_arg(length(rg) == N + 1, "column_solve: r_grid must have N + 1 faces")

        # enclosed mass at the faces, consistent with cs.g = G M(r_face above) / r_face_above^2 (resample_profiles)
        Mf = vcat(cs.g[1] * rg[2]^2 / G - 4π / 3 * cs.ρ[1] * (rg[2]^3 - rg[1]^3), cs.g .* rg[2:end] .^ 2 ./ G)
        
        # gravity at the node (H2 of obliqua_grid_order.md): M(r_n) = M(face below) + (4 pi/3) rho_n (r_n^3 - r_face^3)
        gnode = [G * (Mf[i] + 4π / 3 * cs.ρ[i] * (cs.r_c[i]^3 - rg[i]^3)) / cs.r_c[i]^2 for i in 1:N]
        gcell = g_at === :node ? gnode : cs.g
        
        # non-dimensional node quantities
        rs = cs.r_c ./ R0; ρs = cs.ρ ./ ρ0; gs = gcell ./ g0; μs = cs.μ ./ μ0; Ks = cs.κs ./ μ0; ωs = cs.σ / ω0
        rgs = rg ./ R0
        face = bc_at === :face
        ρₗs = cs.ρl ./ ρ0; Kls = cs.κl ./ μ0; Kds = cs.κd ./ μ0; ηₗs = cs.ηl ./ (μ0 / ω0); ks = cs.k ./ R0^2
        σd = cs.σR ./ ω0
        kw = (G0 = G0, inertial_terms = cs.inertial_terms)
        # layered profile (porosity steps at the input layer faces) or nodal profile (material varies at every node);
        # a column built without column_spec has no such field and counts as layered
        stepped = get(cs, :profile, :layers) === :layers

        # motion matrix of node i in its own slot order (6x6 for solid/fluid, 8x8 for mush)
        function Anode(i)
            ti = t[i]
            # Rayleigh drag (drag = σ_R / ω0) on the tangential motion of the fluid, and with fluid_criterion = :skeleton
            # on the matrix of the solid/mush nodes below it (eq:A_drag)
            (ti === :solid || ti === :fluid) && return get_A(ωs, rs[i], ρs[i], gs[i], μs[i], Ks[i], n; kw..., Y = Y6, drag = σd[i])
            ti === :mush && return get_A(ωs, rs[i], ρs[i], gs[i], μs[i], Ks[i], ρₗs[i], Kls[i], Kds[i], cs.α[i], ηₗs[i],
                                         cs.ϕ[i], ks[i], n; kw..., Y = Y8, drag = σd[i])
            error("column_solve: unknown node type $ti")
        end

        # node slots (6 or 8) of each node, in the order of the state vector y
        sl = [node_slots(ti) for ti in t]
        # degrees of freedom of each node (6 or 8)
        M = length.(sl)                                  # state size of each node (6 or 8)
        # a cell is an interface cell where the material type changes. A porosity change inside a mush is an ordinary
        # cell (W is continuous), but in a layered profile its half-steps still meet at the layer face (atface)
        porojump(i) = t[i] === :mush && t[i+1] === :mush && cs.ϕ[i] != cs.ϕ[i+1]
        isface(i) = t[i] !== t[i+1]
        atface(i) = isface(i) || (stepped && porojump(i))

        # ---- cell equations C_n y_n + D_{n+1} y_{n+1} = 0, split into upper (u) and lower (l) halves. 
        # Block row n of the Henyey system holds the lower half of cell n-1 and the upper half of cell n.
        Cu = Vector{Matrix{ComplexF64}}(undef, N - 1); Du = similar(Cu); Cl = similar(Cu); Dl = similar(Cu)
        faces = Any[]

        # start recursion at the CMB with the first node's motion matrix, then loop over the cells
        A_np = Anode(1)
        for i in 1:N-1
            A_n = A_np; A_np = Anode(i + 1)
            dr = rs[i+1] - rs[i]
            # half-step lengths: dr/2 each in an ordinary cell (trapezoidal rule); in an interface cell the
            # node-to-face distances h- = r_face - r_n and h+ = r_{n+1} - r_face (eq:face_values), so that the face
            # values are extrapolated to the actual material face r_grid[i+1]
            hlo, hhi = (face && atface(i)) ? (rgs[i+1] - rs[i], rs[i+1] - rgs[i+1]) : (0.5 * dr, 0.5 * dr)
            # C_n = I + h- A_n extrapolates y_n up to the face (y- = C_n y_n);
            # D_{n+1} = -I + h+ A_{n+1} extrapolates y_{n+1} down to the face (y+ = -D_{n+1} y_{n+1})
            C = Matrix{ComplexF64}(I, M[i], M[i]) + hlo * A_n
            D = -Matrix{ComplexF64}(I, M[i+1], M[i+1]) + hhi * A_np
            if !isface(i)
                # ordinary cell: every row is continuity (L = G = I); rows 1:M/2 are the upper half (slots U, V, Phi
                # [, P]) and go into block row i, rows M/2+1:M are the lower half (X, Y, Psi [, W]) for block row i+1
                h = M[i] ÷ 2
                Cu[i], Du[i], Cl[i], Dl[i] = C[1:h, :], D[1:h, :], C[h+1:end, :], D[h+1:end, :]
            else
                # material face: the rows are the interface conditions of tab:interface_conditions
                up, low = interface_conditions(t[i], t[i+1]; slip = cs.slip, pores = cs.pores)
                Cu[i], Du[i], Cl[i], Dl[i] = interface_step(C, D, up, low, sl[i], sl[i+1])
                push!(faces, (i = i, lo = t[i], hi = t[i+1], up = up, low = low, C = C, D = D))
            end
        end

        # ---- core rows B_1, applied to y(r_C) = E_1 y_1 with E_1 = I - h_1 A_1, h_1 = r_1 - r_C
        ρc, μc_, κc_ = cs.ρ_core / ρ0, cs.μ_core / μ0, cs.κ_core / μ0

        # CMB radius and gravity at the CMB (core mass only) for bc_at = :face; Obliqua's first node for :node
        # e.g. use core properties for the CMB interface, or the first node's properties for the core as a cell
        rC, gC = face ? (rgs[1], G * Mf[1] / rg[1]^2 / g0) : (rs[1], gs[1])

        # the interface is placed at the bottom of the first node, so the first node's motion matrix is extrapolated down to the CMB
        E1 = face ? Matrix{ComplexF64}(I, M[1], M[1]) - (rs[1] - rgs[1]) * Anode(1) : Matrix{ComplexF64}(I, M[1], M[1])

        # core rows: left null space of the core solutions (a liquid core as inertial where finite)
        B, core_used = core_rows(cs.core, t[1], ωs, rC, ρc, gC, μc_, κc_, n; G0 = G0)
        B = B * E1

        # ---- forward elimination y_i = R_i y_{i+1} (M_n = 0: no sources but the interface terms)
        Rm = Vector{Matrix{ComplexF64}}(undef, N - 1)
        # R_n = -X_n^{-1} Q_n; a direct LU solve, with the minimum-norm pseudo-inverse as fallback for a singular or
        # rectangular X (never reached for the physical columns tested)
        solve(X, Q) = size(X, 1) == size(X, 2) ? (F = lu(X; check = false); issuccess(F) ? -(F \ Q) : -pinv(X) * Q) :
                                                 -pinv(X) * Q
        # block row 1: S_1 = [B_1 E_1; C_1^u], Q_1 = [0; D_2^u]
        Rm[1] = solve([B; Cu[1]], [zeros(ComplexF64, size(B, 1), M[2]); Du[1]])
        for i in 2:N-1
            # block row i: X_i = P_i R_{i-1} + S_i with P_i = [C_{i-1}^l; 0], S_i = [D_i^l; C_i^u], Q_i = [0; D_{i+1}^u]
            Rm[i] = solve([Cl[i-1] * Rm[i-1] + Dl[i-1]; Cu[i]], [zeros(ComplexF64, size(Cl[i-1], 1), M[i+1]); Du[i]])
        end

        # ---- surface rows B_N on y(R) = E_N y_N, E_N = I + h_N A_N, h_N = R - r_N
        Sd = [S[i, i] for i in 1:8]                        # diagonal scale of the standard 8-state order
        std_of = [findfirst(==(s), Y8) for s in 1:8]       # slot s -> standard component
        loc(s) = findfirst(==(s), sl[N])                   # position of slot s in the top node's state

        # the top radius and gravity for the Love numbers, depending on bc_at (extrapolated to the face or taken at the node)
        rtop, gtop = face ? (rg[end], cs.g[end]) : (cs.r_c[end], cs.g[end])
        EN = face ? Matrix{ComplexF64}(I, M[N], M[N]) + (rgs[end] - rs[N]) * Anode(N) : Matrix{ComplexF64}(I, M[N], M[N])
        
        # surface rows B_N and b
        Bt = zeros(ComplexF64, M[N] ÷ 2, M[N]); bt = zeros(ComplexF64, M[N] ÷ 2)

        # row `row` of B_N picks slot s of the physical state: (S z)_s
        put!(row, s) = (Bt[row, loc(s)] = Sd[std_of[s]])
        
        # free surface: X(R) = 0, Y(R) = 0, Psi(R) = (2n+1)/R for a unit degree-n tidal potential
        put!(1, SX); put!(2, SY); put!(3, SΨ); bt[3] = (2n + 1) / rtop
        # porous surface node: closed (impermeable) surface, W(R) = 0
        M[N] == 8 && put!(4, SR)

        # apply extropolation step
        Bt = Bt * EN

        # last block row: X_N = [C_{N-1}^l R_{N-1} + D_N^l; B_N E_N], m_N = [0; b]. Solved directly: the column
        # system has no near-singular surface matrix from a separate segment (the SVD floor of the segment solvers
        # is not needed here, and would bias k).
        X = [Cl[N-1] * Rm[N-1] + Dl[N-1]; Bt]
        z = Vector{Vector{ComplexF64}}(undef, N)
        z[N] = X \ [zeros(ComplexF64, size(Cl[N-1], 1)); bt]
        # back substitution
        for i in N-1:-1:1
            z[i] = Rm[i] * z[i+1]
        end

        # physical profiles in the standard order (U, V, y3, y4, Phi, Psi, P, W); P and W stay 0 at non-porous nodes
        y = zeros(ComplexF64, 8, N)
        for i in 1:N, (j, s) in enumerate(sl[i])
            y[std_of[s], i] = Sd[std_of[s]] * z[i][j]   # renormalize to physical units
        end

        # Love numbers at r_top (R for bc_at = :face): k = Phi(R) - 1 for a unit tidal potential, h = g U(R)
        zR = EN * z[N]      # extrapolate the top node's state to the top radius (face or node)
        kR = Sd[std_of[SΦ]] * zR[loc(SΦ)] - 1
        hR = gtop * Sd[std_of[SU]] * zR[loc(SU)]
        return (k = kR, h = hR, y = y, z = z, faces = faces, slots = sl, types = t, finite = all(isfinite, z[N]),
                r_top = rtop, scale = [Sd[std_of[s]] for s in 1:8], core = core_used, g_node = gnode, z_C = E1 * z[1],
                z_R = zR)
    end


    """
        column_heating(cs, y, n, ω; solid=:displacement, jump=nothing)

    Dissipation of the column solution `y` (column_solve) per solver cell, from `common.get_heating_profile_nodal`: fields
    `shear`, `bulk`, `darcy`, `drag` (shell-mean densities, W m^-3 per |U_T|^2) and the integrated power `P_total`
    (W per |U_T|^2).

    Strain recovery (`common.strain_parts`): fluid nodes always use the displacement-based strain (y4 / mu_f is ill
    conditioned for a fluid); solid and mush nodes use the same (solid = :displacement) or the stress-based node strain
    (solid = :stress). `jump` (default: every change of node type and, in a layered profile, every porosity step
    inside a mush) marks the faces across which no displacement difference is taken.
    """
    function column_heating(cs, y, n::Int, ω::Real; solid::Symbol = :displacement, jump = nothing)
        disp, porous, jmp = heating_masks(cs, solid, jump)
        Eμ, Eκ, El, Ed = get_heating_profile_nodal(y, cs.r_c, cs.r_grid, cs.ρ, cs.g, cs.μ, cs.κs, cs.κd, cs.α, cs.ϕ, cs.κl,
                                                   cs.ηl, cs.k, cs.ρl, cs.σR, porous, n, ω; disp = disp, jump = jmp)
        dvol = 4π / 3 .* (cs.r_grid[2:end] .^ 3 .- cs.r_grid[1:end-1] .^ 3)
        return (shear = Eμ, bulk = Eκ, darcy = El, drag = Ed, solid_method = solid,
                P_total = sum((Eμ .+ Eκ .+ El .+ Ed) .* dvol))
    end


    """
        column_heating_map(cs, y, n, ω, SphericalGrid; solid=:displacement, jump=nothing)

    Heating maps (lat x lon x solver cell) of the column solution, `common.get_heating_map_nodal` with the strain
    recovery of `column_heating`: shear, bulk, Darcy and drag.
    """
    function column_heating_map(cs, y, n::Int, ω::Real, SphericalGrid::NamedTuple; solid::Symbol = :displacement,
                                jump = nothing)
        disp, porous, jmp = heating_masks(cs, solid, jump)
        return get_heating_map_nodal(y, cs.r_c, cs.r_grid, cs.ρ, cs.g, cs.μ, cs.κs, cs.κd, cs.α, cs.ϕ, cs.κl, cs.ηl, cs.k,
                                     cs.ρl, cs.σR, porous, n, ω, SphericalGrid; disp = disp, jump = jmp)
    end

    # strain recovery of the column: displacement-based nodes, porous (8-state) nodes, and the faces across which no
    # displacement difference is taken (default: every change of node type and, in a layered profile, every porosity
    # step inside a mush; in a nodal profile the porosity changes at every node and varies smoothly)
    function heating_masks(cs, solid::Symbol, jump)
        check_option("column_heating", "solid", solid, (:displacement, :stress))
        Nc = length(cs.r_c)
        stepped = get(cs, :profile, :layers) === :layers
        jmp = jump === nothing ? [cs.types[j] !== cs.types[j+1] ||
                                  (stepped && cs.types[j] === :mush && cs.ϕ[j] != cs.ϕ[j+1]) for j in 1:Nc-1] : jump
        return cs.fluid .| (solid == :displacement), cs.types .=== :mush, jmp
    end


    """
        solve_column(cs; n=2, solid=:auto)

    One frequency of the column model: column_solve, heating and energy closure. solid = :auto evaluates both strain
    recoveries for the solid/mush nodes (fluid nodes are always displacement-based) and keeps the one whose closure
    P_total / tidal work is closer to 1 (both are consistent discretisations of the same dissipation; the closure picks
    the one that resolves the boundary layers of this solution better).

    # Returns
    NamedTuple: `k`, `h`, `sol` (column_solve), `heating` (column_heating of the chosen recovery), `P_tide`,
    `closure`, `closure_disp`, `closure_stress`, `passive` (-sign(ω) Im k >= 0).
    """
    function solve_column(cs; n::Int = 2, solid::Symbol = :auto)
        sol = column_solve(cs, n)
        ω = cs.σ
        P_tide = tidal_work(sol.k, n, sol.r_top, ω)
        hD = column_heating(cs, sol.y, n, ω; solid = :displacement)
        passive = sign(ω) * imag(sol.k) <= 0
        if all(cs.fluid) || solid == :displacement
            return (k = sol.k, h = sol.h, sol = sol, heating = hD, P_tide = P_tide, closure = hD.P_total / P_tide,
                    closure_disp = hD.P_total / P_tide, closure_stress = NaN, passive = passive)
        end
        hS = column_heating(cs, sol.y, n, ω; solid = :stress)
        cD, cS = hD.P_total / P_tide, hS.P_total / P_tide
        good(c) = isfinite(c) && c > 0
        h = solid == :stress ? hS :
            !good(cS) ? hD : !good(cD) ? hS : abs(log(cD)) <= abs(log(cS)) ? hD : hS
        return (k = sol.k, h = sol.h, sol = sol, heating = h, P_tide = P_tide, closure = h.P_total / P_tide,
                closure_disp = cD, closure_stress = cS, passive = passive)
    end

end

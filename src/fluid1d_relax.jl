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
#            in the tangential momentum row. The fluid layers are those from the surface down below a viscosity
#            cutoff (fluid_mask), the same at every forcing frequency.
#
# Interfaces (sec:general_interface): every face is a list of linear conditions L y- = G y+ on the face values
# y- = C_n y_n and y+ = -D_{n+1} y_{n+1} (eq:face_values, eq:interface_LG), which turns the cell equation into
# (L C_n) y_n + (G D_{n+1}) y_{n+1} = 0 (eq:interface_cell). Paper subset (tab:interface_conditions):
#   free slip at every solid/mush <-> fluid face and at a liquid core (Y = 0 on the side(s) that carry shear, V free),
#   closed pores at solid <-> mush (R = 0 on the porous side, P free),
#   open pores at mush <-> fluid (U_f = U_s + phi R, P + X_f = 0),
#   melt flux continuity phi- R- = phi+ R+ where the porosity changes between two mush cells (everything else continuous).
# Core: a static liquid core is the 4x2 solution (eq:IC_liquid) and the CMB is an interface like any other; above
# omega_sw the "inertial-liquid" core of Obliqua is used instead (its 6x3 I_C has the same free-slip structure).
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
#   y = (y1, ..., y8) = (U, V, sigma_rr, sigma_rtheta, Phi, Psi, P, R)
#   U, V   radial and tangential displacement [m]
#   X = y3 = sigma_rr, Y = y4 = sigma_rtheta [Pa]
#   Phi = y5 potential perturbation, Psi = y6 its gradient term (Obliqua's sign convention, k = y5(R) - 1)
#   P = y7 pore pressure [Pa], R = y8 radial displacement of the melt relative to the matrix [m]
# The solver works in the "slot" order of the appendix (eq:y_order), U V Phi P | X Y Psi R, so that the upper half
# (slots 1-4) holds the quantities bounded from below and the lower half (slots 5-8) the tractions:
#   storage slot   1  2  3    4  5  6  7    8
#   component      U  V  Phi  P  X  Y  Psi  R
# A solid or fluid node carries only the six slots [1, 2, 3, 5, 6, 7]; a mush node all eight. The permutation is
# passed to Obliqua's get_A as Y (Y[k] = slot of standard component k): Y6 = [1, 2, 4, 5, 3, 6] for the 6-state
# order (U, V, Phi, X, Y, Psi) and Y8 = [1, 2, 5, 6, 3, 7, 4, 8] for the 8-state order.

module fluid1d_relax

    include("constants.jl")
    include("common.jl")
    using .constants
    using .common

    using LinearAlgebra
    using Printf

    export column_spec, column_solve, column_heating, column_heating_map, solve_column, face_residuals, tidal_work
    export interface_conditions, interface_step, liquid_core_rows, A_drag, node_slots, refine_layers, fluid_mask
    export boundary_layer_resolution

    # storage slots: U 1, V 2, Phi 3, P 4, X 5, Y 6, Psi 7, R 8
    const SU, SV, SΦ, SP, SX, SY, SΨ, SR = 1, 2, 3, 4, 5, 6, 7, 8

    # get_A permutations: Y[k] is the slot that standard component k (U, V, y3, y4, y5, y6[, P, R]) occupies
    const Y6 = [1, 2, 4, 5, 3, 6]          # 6-state order U, V, Phi, X, Y, Psi
    const Y8 = [1, 2, 5, 6, 3, 7, 4, 8]    # 8-state order U, V, Phi, P, X, Y, Psi, R

    """
        node_slots(t)

    Storage slots carried by a node of type `t`: all eight for a mush, the six non-porous slots [1, 2, 3, 5, 6, 7]
    (U, V, Phi, X, Y, Psi) for a solid or fluid node. Each node keeps its own size in the recursion.
    """
    node_slots(t::Symbol) = t === :mush ? collect(1:8) : [SU, SV, SΦ, SX, SY, SΨ]

    """Unit vector e_s of C^8 (8-slot storage order)."""
    unit8(s::Int) = (v = zeros(8); v[s] = 1.0; v)


    # ------------------------------------------------------------------------------------------------ motion matrix

    """
        A_drag(ω, r, ρ, g, μ, K, n, σd; G0, inertial_terms, Y)

    6x6 motion matrix of a fluid node with absolute Rayleigh drag, non-dimensional:

        A_f = A(mu_f, K_f) + i w rho sigma_R e_Y e_V^T.

    Obliqua's get_A puts -rho w^2 in the (Y, V) entry, i.e. the tangential momentum row dY/dr has the term
    -rho w^2 V from the tangential inertia. With drag the inertia becomes rho w (w - i sigma_R), so the entry gains
    +i w rho sigma_R. Only the tangential row changes: the radial inertia (the (X, U) entry) is left as it is, so the
    thin-layer limit is the drag of the Laplace tidal equations. With sigma_R = 0 this is exactly Obliqua's matrix.
    The drag does work -1/2 rho sigma_R w^2 integral |u_h|^2 dV <= 0, so the response is passive for sigma_R >= 0.

    # Arguments
    - `ω, r, ρ, g, μ, K`  : Non-dimensional frequency, radius, density, gravity and complex moduli of the node.
    - `n::Int`            : Spherical-harmonic degree.
    - `σd`                : Non-dimensional drag rate sigma_R / omega_0 (0 for solid nodes).
    """
    function A_drag(ω, r, ρ, g, μ, K, n, σd; G0, inertial_terms, Y)
        A = get_A(ω, r, ρ, g, μ, K, n; G0 = G0, inertial_terms = inertial_terms, Y = Y)
        # row of Y (= y4, tangential stress), column of V (= y2): tangential inertia -rho w^2 -> -rho w (w - i sigma_R)
        iszero(σd) || (A[Y[4], Y[2]] += 1im * ω * ρ * σd)
        return A
    end


    # ------------------------------------------------------------------------------------------- sub-layer refinement

    """
        refine_layers(n_sub, r, ρ, η, μc, κc, κdc, κl, ϕ, k; κu=nothing) -> NamedTuple

    Split every input layer into `n_sub` equal radial sub-layers and interpolate the layer values linearly in radius
    between the layer centres (flat beyond the first and last centre). The fluid mask (`fluid_mask`) is evaluated
    on the sub-layers, so an interface moves with the layer values instead of jumping by a whole layer when one layer
    changes class.

    Interpolation: density and porosity linearly; viscosity, melt and undrained bulk moduli and permeability in log (linearly where
    an end value is zero); the real and imaginary parts of the complex moduli separately, in log of the magnitude where
    both end values share a sign and linearly otherwise, so passive moduli stay passive and moduli at -ω stay the
    conjugates of those at ω. The Biot coefficient is recomputed from the interpolated moduli by the caller.

    # Returns
    NamedTuple (r, ρ, η, μc, κc, κdc, κl, ϕ, k, κu, parent) on the sub-layers (`κu` nothing when not given); `parent` is
    the input layer of each.
    """
    function refine_layers(n_sub::Integer, r::AbstractVector, ρ::AbstractVector, η::AbstractVector,
                           μc::AbstractVector, κc::AbstractVector, κdc::AbstractVector, κl::AbstractVector,
                           ϕ::AbstractVector, k::AbstractVector; κu = nothing)
        # check arguments
        n_sub >= 1 || throw(ArgumentError("refine_layers: n_sub = $n_sub must be >= 1"))
        
        # sub-layer centres and weights for interpolation
        Nl = length(ρ)
        rc0 = 0.5 .* (r[1:end-1] .+ r[2:end])
        rs = vcat([collect(range(r[i], r[i+1]; length = n_sub + 1))[1:end-1] for i in 1:Nl]..., [r[end]])
        rcs = 0.5 .* (rs[1:end-1] .+ rs[2:end])
        # input layer of each sub-layer (1:Nl) for interpolation
        parent = repeat(1:Nl, inner = n_sub)

        # bracketing centres (a, b) and weight w of each sub-layer centre; flat outside [rc0[1], rc0[end]]
        a = similar(parent); b = similar(parent); w = zeros(Float64, length(rcs))

        # find the bracketing layers for each sub-layer centre
        for (s, x) in enumerate(rcs)
            j = clamp(searchsortedlast(rc0, x), 1, Nl)
            if x <= rc0[1] || Nl == 1
                a[s] = b[s] = 1
            elseif x >= rc0[end]
                a[s] = b[s] = Nl
            else
                a[s] = j; b[s] = j + 1; w[s] = (x - rc0[j]) / (rc0[j+1] - rc0[j])
            end
        end

        # linear interpolation of the input values to the sub-layer centres, with flat extrapolation beyond the first and last centre
        lin(v) = [(1 - w[s]) * v[a[s]] + w[s] * v[b[s]] for s in eachindex(rcs)]
        # log interpolation where both ends are positive, linear otherwise (a zero stays a zero at its centre)
        lg(x, y, t) = t == 0 ? x : t == 1 ? y : (x > 0 && y > 0) ? exp((1 - t) * log(x) + t * log(y)) : (1 - t) * x + t * y
        lgv(v) = [lg(v[a[s]], v[b[s]], w[s]) for s in eachindex(rcs)]
        # complex parts in log of the magnitude where both ends share a sign, so the moduli at -ω stay the conjugates
        lgs(x, y, t) = (x < 0 && y < 0) ? -lg(-x, -y, t) : lg(x, y, t)
        lgc(v) = [complex(lgs(real(v[a[s]]), real(v[b[s]]), w[s]), lgs(imag(v[a[s]]), imag(v[b[s]]), w[s]))
                  for s in eachindex(rcs)]
        return (r = rs, ρ = lin(ρ), η = lgv(η), μc = lgc(μc), κc = lgc(κc), κdc = lgc(κdc), κl = lgv(κl), ϕ = lin(ϕ),
                k = lgv(k), κu = κu === nothing ? nothing : lgv(κu), parent = parent)
    end


    # ------------------------------------------------------------------------------------------------ column spec

    """
        column_spec(ω, r, ρ, η, μc, κc, κdc, αB, κl, ηl, ϕ, k, m_core, ρ_core, μ_core, κ_core; kwargs...)

    Frequency-dependent column of the whole interior (core to surface) on Obliqua's resampled grid, as a NamedTuple
    for `column_solve` and `column_heating`.

    The layers from the surface down below the viscosity `visc_fluid` are fluid (`fluid_mask`): they keep their shear
    modulus, take the undrained (Gassmann) modulus of their frame, grain and melt moduli as their only bulk modulus (the
    mush limit without Darcy flow; the melt modulus for phi = 1), no porosity, and the drag sigma_R. Porosity below                 ### Check overwritten values
    `porosity_thresh` is set to zero (as in the mush solver), so a layer with phi < porosity_thresh is solid. With `κu`
    (the undrained modulus of the interior rheology, per layer) that modulus replaces the Gassmann modulus built here
    wherever it is positive; the interior computes it with the melt volume fraction.

    # Arguments (per input layer unless noted; r has N_l + 1 boundaries)
    - `ω`               : Forcing frequency [rad s^-1], nonzero; the moduli are those at ω (conjugate for ω < 0).
    - `r, ρ, η`         : Layer boundaries [m], density [kg m^-3], shear viscosity [Pa s].
    - `μc, κc, κdc, αB` : Complex shear, grain bulk and drained bulk moduli [Pa] at ω, and the Biot coefficient.
    - `κl`              : Melt bulk modulus [Pa] per layer.
    - `ηl`              : Melt viscosity [Pa s] (scalar), used in the Darcy terms of the mush.
    - `ϕ, k`            : Porosity [-] and permeability [m^2].
    - `m_core, ρ_core`  : Core mass [kg] and density [kg m^-3].
    - `μ_core, κ_core`  : Core moduli [Pa] (only used by solid/inertial core solutions).

    # Keyword Arguments
    - `σ_R=1e-4`             : Rayleigh drag rate in the fluid layers [s^-1], >= 0.
    - `dr_min, dr_max`       : Grid spacing at the surface and at the CMB [m] (Obliqua's exponential grid).
    - `core="liquid"`        : Configured core type.
    - `ω_sw=1e-5`            : A "liquid" core becomes "inertial-liquid" for |ω| >= ω_sw [rad s^-1].
    - `inertial_terms=true`  : Inertia in the motion matrices.
    - `porosity_thresh=0.0`  : Porosity below which a layer is non-porous.
    - `slip=:freeslip`       : :freeslip or :noslip at solid/mush <-> fluid faces.
    - `pores=:open`          : :open or :closed at mush <-> fluid faces (solid <-> mush is always closed).
    - `visc_fluid=1e8`       : Viscosity [Pa s] below which the layers connected to the surface are fluid, at every
                               frequency. Layers with η >~ 1e10 Pa s as fluid break the energy closure of the heating.
    - `n_sub=1`              : Interpolated sub-layers per input layer (`refine_layers`).
    - `κu=nothing`           : Undrained bulk modulus per input layer [Pa] (entries <= 0: Gassmann of the layer here).

    A layer below the porosity threshold is solved as non-porous with the undrained (Gassmann) bulk modulus of its
    melt-bearing frame, K_u = K_d + alpha^2 / S (common.gassmann), as both its grain and drained modulus. A poorly
    connected mush with little melt is undrained at tidal frequencies, so the bulk stiffness is continuous where a
    layer crosses the threshold; without melt K_u = K_s.
    """
    function column_spec(ω::Real, r::AbstractVector, ρ::AbstractVector, η::AbstractVector, μc::AbstractVector,
                         κc::AbstractVector, κdc::AbstractVector, αB::AbstractVector, κl::AbstractVector, ηl::Real,
                         ϕ::AbstractVector, k::AbstractVector, m_core::Real, ρ_core::Real, μ_core::Number,
                         κ_core::Number; σ_R::Real = 1e-4, dr_min::Integer = 300, dr_max::Integer = 3000,
                         core::String = "liquid", ω_sw::Real = 1e-5, inertial_terms::Bool = true,
                         porosity_thresh::Real = 0.0, slip::Symbol = :freeslip, pores::Symbol = :open,
                         visc_fluid::Real = 1e8, n_sub::Integer = 1, κu = nothing)
        # check arguments
        ω != 0 || throw(ArgumentError("column_spec: ω must be nonzero (evaluate the relaxed limit at a small ω)"))
        visc_fluid > 0 || throw(ArgumentError("column_spec: visc_fluid = $visc_fluid must be > 0"))
        σ_R >= 0 || throw(ArgumentError("column_spec: σ_R = $σ_R must be >= 0 (passivity)"))
        slip in (:freeslip, :noslip) || throw(ArgumentError("column_spec: slip must be :freeslip or :noslip"))
        pores in (:open, :closed) || throw(ArgumentError("column_spec: pores must be :open or :closed"))
        Nl = length(ρ)
        length(r) == Nl + 1 || throw(ArgumentError("column_spec: r needs length(ρ) + 1 boundaries"))
        # the exponential grid of resample_profiles needs ln(dr_max/dr_min) > 0
        dr_max > dr_min > 0 || throw(ArgumentError("column_spec: need dr_max > dr_min > 0, got $dr_min, $dr_max"))
        n_sub >= 1 || throw(ArgumentError("column_spec: n_sub = $n_sub must be >= 1"))
        
        # refine the layers around the interfaces if applicable
        if n_sub > 1
            # sub-layers with interpolated values, so class boundaries move continuously (refine_layers)
            sub = refine_layers(n_sub, r, ρ, η, μc, κc, κdc, κl, ϕ, k; κu = κu)
            r, ρ, η, μc, κc, κdc, κl, ϕ, k, κu = sub.r, sub.ρ, sub.η, sub.μc, sub.κc, sub.κdc, sub.κl, sub.ϕ, sub.k, sub.κu
            αB = 1 .- κdc ./ κc
            Nl = length(ρ)
        end

        # the low-viscosity layers connected to the surface are fluid, the same at every frequency
        fluid = fluid_mask(η; visc_fluid = visc_fluid)
        any(i -> fluid[i] && iszero(μc[i]), 1:Nl) &&
            throw(ArgumentError("column_spec: a fluid layer has a zero shear modulus (give the melt μ and η)"))
        
        # resample the input profiles to Obliqua's exponential grid
        rr = Vector{prec}(r)
        ρa = Vector{prec}(ρ); ηa = Vector{prec}(η); ϕa = Vector{prec}(ϕ); ka = Vector{prec}(k)
        μa = Vector{precc}(μc); κa = Vector{precc}(κc); κda = Vector{precc}(κdc); αa = Vector{precc}(αB)

        # apply the porosity threshold, as in run_solid1d_mush_relax. A layer below it is solved as non-porous with the
        # undrained (Gassmann) modulus of its frame as its only bulk modulus (K_u = K_s without melt), alpha = 0
        κu === nothing || length(κu) == Nl || throw(ArgumentError("column_spec: κu needs one entry per layer"))
        # undrained modulus of layer i: the interior's where given, else Gassmann with ϕ as porosity
        undrained(i) = (κu !== nothing && κu[i] > 0) ? precc(κu[i]) : gassmann(κda[i], κa[i], prec(κl[i]), ϕa[i])
        nopore = ϕa .< porosity_thresh
        for i in findall(nopore)
            Ku = undrained(i)
            κa[i] = Ku; κda[i] = Ku; αa[i] = 0.0
        end
        ϕa[nopore] .= 0.0

        # fluid layers: undrained bulk modulus of the frame (Gassmann), no porosity, Rayleigh drag on the tangential
        for i in findall(fluid)
            # undrained bulk modulus of the layer (kappa_l at phi = 1) as grain and drained modulus; no porosity                    # Check this 
            κa[i] = undrained(i)
            κda[i] = κa[i]; αa[i] = 0.0
            ϕa[i] = 0.0; ka[i] = 0.0
        end

        # melt bulk modulus and viscosity per layer
        κla = Vector{prec}(vcat(κl, κl[end]))
        ηla = fill(prec(ηl), Nl + 1)

        # Obliqua's exponential grid: cell faces r_grid, material constant per cell, g_grid = G M(face above) / r^2 at the outer faces
        r_grid, ρg, ηg, μg, κsg, κlg, κdg, αg, ηlg, ϕg, kg, g_grid, M_tot =
            resample_profiles(rr, ρa, ηa, μa, κa, κla, κda, αa, ηla, ϕa, ka, prec(m_core), Int(dr_min), Int(dr_max))
        
        # cell centers of the resampled grid
        r_c = 0.5 .* (r_grid[1:end-1] .+ r_grid[2:end])

        # node -> input layer (resample_profiles' own rule), and the fluid mask of the nodes
        lay = [clamp(searchsortedfirst(rr, x) - 1, 1, Nl) for x in r_c]
        fl = [Bool(fluid[L]) for L in lay]

        # Rayleigh drag profile for fluid layers (uniform)
        σR = [fl[q] ? Float64(σ_R) : 0.0 for q in eachindex(r_c)]

        # node types: fluid, mush or solid
        types = [fl[q] ? :fluid : (ϕg[q] > 0 && kg[q] > 0) ? :mush : :solid for q in eachindex(r_c)]

        # core switch: Obliqua's "inertial-liquid" core returns NaN below ~5e-6 rad/s in Float64; the static
        # "liquid" core is used below omega_sw
        core_type = (core == "liquid" && abs(ω) >= ω_sw) ? "inertial-liquid" : core
        return (σ = Float64(ω), r_grid = r_grid, r_c = r_c, ρ = ρg, g = g_grid, μ = μg, κs = κsg, κl = κlg, κd = κdg,
                α = αg, ηl = ηlg, ϕ = ϕg, k = kg, ρl = ρg .* ϕg, scales = [rr[end], M_tot, G], M_tot = M_tot,
                ρ_core = Float64(ρ_core), μ_core = ComplexF64(μ_core), κ_core = ComplexF64(κ_core), core = core_type,
                inertial_terms = inertial_terms, fluid = fl, σR = σR, types = types, lay = lay, slip = slip,
                pores = pores, fluid_layers = fluid)
    end


    """
        fluid_mask(η; visc_fluid=1e8) -> Vector{Bool}

    Fluid mask of the layers (ordered bottom to top): the contiguous layers from the surface down whose viscosity η is
    below `visc_fluid` [Pa s]. They get the fluid treatment (the undrained bulk modulus, no melt variables, the Rayleigh
    drag) and keep their own shear modulus. The mask does not depend on the forcing frequency, so every tidal mode of a
    snapshot sees the same layering. For the tidal frequencies and the melt and mush shear moduli of interest a layer
    below about 1e8 Pa s also relaxes faster than the tide (Maxwell time eta / mu below 1 / omega).

    Only the ocean connected to the surface is fluid. A buried layer below the cutoff stays a solid or mush layer: with a
    non-monotonic viscosity near the cutoff it would otherwise form a second fluid region under a thin mush sliver, which
    changes k by tens of per cent at every frequency.
    """
    function fluid_mask(η::AbstractVector; visc_fluid::Real = 1e8)::Vector{Bool}
        out = falses(length(η))
        i = length(η)
        while i >= 1 && η[i] < visc_fluid
            out[i] = true
            i -= 1
        end
        return collect(Bool, out)
    end


    # ------------------------------------------------------------------------------------------------ interfaces

    """
        interface_conditions(lo, hi; slip=:freeslip, pores=:open, ϕlo=0, ϕhi=0) -> (up, low)

    Conditions l . y- = g . y+ (pairs (l, g) of 8-slot vectors) at the face between a node of type `lo` (below) and
    `hi` (above). A condition keeps the row of the component whose continuity it replaces:

    - `up`  : the M-/2 rows of block row n, i.e. the upper-half slots of the node below: U, V, Phi and, for a mush
              below, slot 4 (P).
    - `low` : the M+/2 rows of block row n+1, i.e. the lower-half slots of the node above: X, Y, Psi and, for a mush
              above, slot 8 (R).

    In total x = (M- + M+)/2 conditions. A one-sided porous pair (P, R) gives one condition, placed in the row of its
    free member: at the top of a mush (slot 4 of block row n) and at its base (slot 8 of block row n+1). `ϕlo`/`ϕhi` 
    are the porosities of the porous side for the open-pore melt flux.

    Mush | mush with a porosity change (lo === hi === :mush): R is the displacement of the melt relative to the matrix,                     # Is this because we pass porosity phi and not melt mass fraction v?
    and the fluid-content equation (row R of the 8x8 matrix, zeta = -div(phi w)) conserves the volume flux phi R. At a
    porosity jump the melt flux is continuous, phi- R- = phi+ R+, while R itself jumps; this condition takes the R row
    (slot 8 of block row n+1). All other components are continuous: U, V (welded matrix), Phi, the total tractions X, Y,
    the pore pressure P (connected pores, finite Darcy flux), and Psi, which contains the full mass flux
    4 pi G (rho U + rho_l phi R) and is therefore continuous at any density or porosity jump. Without this face the
    recursion imposes R- = R+, which creates or absorbs melt at every porosity jump and breaks the energy closure in
    proportion to the melt flow (Darcy flux), independently of the grid spacing.
    """
    function interface_conditions(lo::Symbol, hi::Symbol; slip::Symbol = :freeslip, pores::Symbol = :open,
                                  ϕlo::Real = 0.0, ϕhi::Real = 0.0)
        z8 = zeros(8)
        # continuity of slot s: e_s . y- = e_s . y+ (the row of an ordinary cell, L = G = I)
        cont(s) = (unit8(s), unit8(s))
        # initialize the rows of the continuity conditions
        rows = Dict(s => cont(s) for s in (SU, SV, SΦ, SX, SY, SΨ))
        # the 4th and 8th rows of block rows n and n+1 are only present for a mush below or above, respectively
        r4 = r8 = nothing

        # helper to check if a node type is fluid
        isfl(t) = t === :fluid

        # check the interface type
        if lo === hi
            # mush | mush with a porosity change: P continuous, melt flux phi- R- = phi+ R+ in place of R continuity.                           # again maybe that can be removed once we pass melt mass fraction instead of porosity, check Hay 2025, the original source
            # (Same-material cells without a porosity change are ordinary cells in column_solve and never get here.)
            lo === :mush && (r4 = cont(SP); r8 = (ϕlo * unit8(SR), ϕhi * unit8(SR)))
        else
            # catch fluid | fluid interfaces, which are not allowed
            isfl(lo) && isfl(hi) && error("interface_conditions: two fluid types do not meet")

            # interface: mush | ...
            if lo === :mush                                               # top of a porous layer (slot 4, block row n)
                # ... | fluid (open pores)
                if pores === :open && isfl(hi)
                    # open pores: pore pressure = fluid pressure, P- + X+ = 0 (P- = -X+, the fluid's sigma_rr)
                    r4 = (unit8(SP), -unit8(SX))
                    # the melt leaving the pores enters the fluid: U- + phi R- = U+ (replaces U continuity)
                    rows[SU] = (unit8(SU) + ϕlo * unit8(SR), unit8(SU))
                # ... | solid/fluid (closed pores)
                else
                    # closed pores: no Darcy flux through the face, R- = 0 (P- is free)
                    r4 = (unit8(SR), z8)
                end
            end
            # interface: ... | mush
            if hi === :mush                                               # base of a porous layer (slot 8, block row n+1)
                # fluid | ... (open pores)
                if pores === :open && isfl(lo)
                    # open pores, mirror image: P+ + X- = 0 and U- = U+ + phi R+
                    r8 = (-unit8(SX), unit8(SP))
                    rows[SU] = (unit8(SU), unit8(SU) + ϕhi * unit8(SR))
                # solid/fluid | ... (closed pores)
                else
                    # closed pores: R+ = 0 (P+ is free)
                    r8 = (z8, unit8(SR))
                end
            end
            # interface: solid/mush/fluid | solid/mush/fluid (free slip) (symmetric), no fluid | fluid interfaces
            if slip === :freeslip && isfl(lo) != isfl(hi)
                # free slip at a solid/mush <-> fluid face: zero shear traction on both sides, V may jump.
                # Y- = 0 replaces V continuity (V row of block row n) ...
                rows[SV] = (unit8(SY), z8)
                # ... and Y+ = 0 replaces Y continuity (Y row of block row n+1)
                rows[SY] = (z8, unit8(SY))
            end
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
        liquid_core_rows(Ic6, t1) -> (B, K, J, Yc)

    The CMB of a static liquid core as an interface ("Fluid core -> solid, mush or fluid (free slip)").                         # Is this not identical to just calling nullspace on the 6x3 I_C? 
    The core is the 4-state (U, X, Phi, Psi) with its two regular solutions Yc (4x2,
    eq:IC_liquid): columns 1 and 3 of Obliqua's liquid I_C (get_Ic in the 6-slot order U, V, Phi, X, Y, Psi),
    restricted to those components; column 2 of that I_C is the pure tangential slip (V only).

    Conditions on the mantle-side face value y+ = y(r_C) (slots of node type t1) and the core weights a:
        S_c y+ - Yc a = 0     continuity of U, X, Phi, Psi (4 rows)
        e_Y . y+ = 0          free slip: no shear traction on the mantle side (V free)
        e_R . y+ = 0          mush only: closed pores at the CMB (P free)
    written as K y+ + J a = 0 with J (m x 2). Eliminating a with N, a basis of the left null space of J (N J = 0),
    gives B_1 = N K with m - 2 = M/2 rows (3 for a solid/fluid node, 4 for a mush).
    """
    function liquid_core_rows(Ic6, t1::Symbol)
        # rows U (slot 1), X (slot 4), Phi (slot 3), Psi (slot 6) of the 6-slot I_C; columns 1 and 3 (regular solutions)
        Yc = Ic6[[1, 4, 3, 6], [1, 3]]
        sl = node_slots(t1)
        loc(s) = findfirst(==(s), sl)
        rowsK = Vector{Vector{ComplexF64}}(); rowsJ = Vector{Vector{ComplexF64}}()
        for (k, s) in enumerate((SU, SX, SΦ, SΨ))
            # mantle component s minus the core's: y+_s - (Yc a)_k = 0
            v = zeros(ComplexF64, length(sl)); v[loc(s)] = 1
            push!(rowsK, v); push!(rowsJ, -Yc[k, :])
        end
        # free slip (Y = 0) and, for a mush, closed pores (R = 0): conditions on the mantle side only
        extra = t1 === :mush ? (SY, SR) : (SY,)
        for s in extra
            v = zeros(ComplexF64, length(sl)); v[loc(s)] = 1
            push!(rowsK, v); push!(rowsJ, zeros(ComplexF64, 2))
        end
        K = permutedims(reduce(hcat, rowsK)); J = permutedims(reduce(hcat, rowsJ))
        # rows of N span the left null space of J, so N (K y + J a) = N K y = 0 for every core state
        Nl = permutedims(nullspace(permutedims(J)))
        return (B = Nl * K, K = K, J = J, Yc = Yc)
    end


    # ------------------------------------------------------------------------------------------------ column solve

    """
        column_solve(cs, n=2; core_route=:auto, bc_at=:face, g_at=:node)

    Response of the column `cs` (column_spec) to a unit degree-n tidal potential, with every material face closed by
    interface_conditions/interface_step. One Henyey recursion from the CMB to the surface.

    # Keyword Arguments
    - `core_route=:auto` : :interface builds B_1 from the CMB as an interface with the 4x2 liquid core
                           (liquid_core_rows; static "liquid" core only); :ic takes B_1 as the left null space of
                           Obliqua's 6x3 I_C (eq:B1_nullspace), with the free pore pressure column added for a mush on
                           the core. :auto uses :interface for a static liquid core and :ic otherwise ("inertial-liquid",
                           "solid", "inertial"). For free slip the two routes span the same rows.
    - `bc_at=:face`      : CMB rows on y(r_C) = (I - h_1 A_1) y_1, surface rows on y(R) = (I + h_N A_N) y_N, k read at
                           R, material faces at their grid radius (second order). :node is Obliqua's convention
                           (boundary rows at the end nodes, faces at node midpoints; first order).
    - `g_at=:node`       : gravity at each node from the enclosed mass (second order); :outer is Obliqua's g of the
                           outer cell face (first order).

    # Returns
    NamedTuple with `k`, `h` (Love numbers at r_top), `y` (8 x N physical profiles, standard order U, V, y3, y4,
    Phi, Psi, P, R), `z` (non-dimensional node states), `faces` (for face_residuals), `slots`, `types`, `finite`,
    `r_top`, `scale`, `core` (core weights and state for the interface route), `g_node`, `z_C`, `z_R`.
    """
    function column_solve(cs, n::Int = 2; core_route::Symbol = :auto, bc_at::Symbol = :face, g_at::Symbol = :node)
        # check arguments
        bc_at in (:face, :node) || error("column_solve: bc_at must be :face or :node")
        g_at in (:node, :outer) || error("column_solve: g_at must be :node or :outer")

        # non-dimensionalization and node quantities
        t = cs.types
        N = length(t)
        N >= 2 || error("column_solve: the column needs at least two nodes")

        # non-dimensionalization (common.get_scales): R0 = planet radius, M0 = planet mass, G0 = G;
        # S is the diagonal state scale, y_physical = S z (8-state order, used for 6-state nodes through their slots)
        R0, M0, ω0, ρ0, G0, g0, μ0, S, Sinv = get_scales(cs.scales[1], cs.scales[2], cs.scales[3]; Y = Y8)
        
        # check grid
        rg = cs.r_grid
        length(rg) == N + 1 || error("column_solve: r_grid must have N + 1 faces")

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

        # motion matrix of node i in its own slot order (6x6 for solid/fluid, 8x8 for mush)
        function Anode(i)
            ti = t[i]
            ti === :solid && return get_A(ωs, rs[i], ρs[i], gs[i], μs[i], Ks[i], n; kw..., Y = Y6)
            ti === :fluid && return A_drag(ωs, rs[i], ρs[i], gs[i], μs[i], Ks[i], n, σd[i]; kw..., Y = Y6)
            ti === :mush && return get_A(ωs, rs[i], ρs[i], gs[i], μs[i], Ks[i], ρₗs[i], Kls[i], Kds[i], cs.α[i], ηₗs[i],
                                         cs.ϕ[i], ks[i], n; kw..., Y = Y8)
            error("column_solve: unknown node type $ti")
        end

        # node slots (6 or 8) of each node, in the order of the state vector y
        sl = [node_slots(ti) for ti in t]
        # degrees of freedom of each node (6 or 8)
        M = length.(sl)                                  # state size of each node (6 or 8)
        # a cell is an interface cell where the material type changes, or where the porosity changes inside a mush
        isface(i) = t[i] !== t[i+1] || (t[i] === :mush && cs.ϕ[i] != cs.ϕ[i+1])

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
            hlo, hhi = (face && isface(i)) ? (rgs[i+1] - rs[i], rs[i+1] - rgs[i+1]) : (0.5 * dr, 0.5 * dr)
            # C_n = I + h- A_n extrapolates y_n up to the face (y- = C_n y_n);
            # D_{n+1} = -I + h+ A_{n+1} extrapolates y_{n+1} down to the face (y+ = -D_{n+1} y_{n+1})
            C = Matrix{ComplexF64}(I, M[i], M[i]) + hlo * A_n
            D = -Matrix{ComplexF64}(I, M[i+1], M[i+1]) + hhi * A_np
            if !isface(i)
                # ordinary cell: every row is continuity (L = G = I); rows 1:M/2 are the upper half (slots U, V, Phi
                # [, P]) and go into block row i, rows M/2+1:M are the lower half (X, Y, Psi [, R]) for block row i+1
                h = M[i] ÷ 2
                Cu[i], Du[i], Cl[i], Dl[i] = C[1:h, :], D[1:h, :], C[h+1:end, :], D[h+1:end, :]
            else
                # material face (or porosity jump in a mush): the rows are the interface conditions of
                # tab:interface_conditions, plus melt flux continuity at a porosity jump
                up, low = interface_conditions(t[i], t[i+1]; slip = cs.slip, pores = cs.pores, ϕlo = cs.ϕ[i],
                                               ϕhi = cs.ϕ[i+1])
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

        # the core route: :interface for a static liquid core, :ic for an inertial-liquid or solid core, :auto for
        # automatic selection
        route = core_route === :auto ? (cs.core == "liquid" ? :interface : :ic) : core_route

        # check if liquid core
        liquid = cs.core in ("liquid", "inertial-liquid")
        core_info = nothing
        # if a liquid core, apply the fluid interface before returning
        if cs.core == "liquid"
            Ic6 = get_Ic(ωs, rC, ρc, gC, μc_, κc_, cs.core, n; G0 = G0, Y = Y6)
            core_info = liquid_core_rows(Ic6, t[1])
        end
        # check if liquid core if an interface (free slip) is needed
        if route === :interface
            core_info === nothing && error("column_solve: core_route = :interface needs a static liquid core")
            # the CMB as an interface with the 4x2 core solution (liquid_core_rows)
            B = core_info.B
        # all other cores don't need an interface condition
        else
            # left null space of I_C (eq:B1_nullspace): y(r_C) = I_C c for some c  <=>  B_1 y(r_C) = 0
            Ic6 = get_Ic(ωs, rC, ρc, gC, μc_, κc_, cs.core, n; G0 = G0, Y = Y6)
            if M[1] == 6
                Ic = Ic6
            else
                # mush on the core: embed I_C in the 8 slots and add the free pore pressure column e_P; the
                # columns have R = 0, so this imposes zero Darcy displacement at the CMB with P free
                Ic = zeros(ComplexF64, 8, 3); Ic[[SU, SV, SΦ, SX, SY, SΨ], :] .= Ic6
                Ic = hcat(Ic, unit8(SP))
            end
            B = permutedims(nullspace(transpose(Ic)))
        end
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
        # porous surface node: closed (impermeable) surface, R(R) = 0
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

        # physical profiles in the standard order (U, V, y3, y4, Phi, Psi, P, R); P and R stay 0 at non-porous nodes
        y = zeros(ComplexF64, 8, N)
        for i in 1:N, (j, s) in enumerate(sl[i])
            y[std_of[s], i] = Sd[std_of[s]] * z[i][j]   # renormalize to physical units
        end

        # infer the core state and weights a from the CMB interface, if any; the misfit is the residual of the core's
        # least-squares problem K y(r_C) + J a = 0, divided by the maximum of |K y(r_C)|
        core = nothing
        if core_info !== nothing
            # core weights a from K y(r_C) + J a = 0 (least squares; exact when the rows hold), and the core state
            zC = E1 * z[1]
            a = -(pinv(core_info.J) * (core_info.K * zC))
            state = core_info.Yc * a                         # (U, X, Phi, Psi) of the core at the CMB, non-dimensional
            core = (a = a, state = state .* [Sd[std_of[s]] for s in (SU, SX, SΦ, SΨ)],
                    misfit = norm(core_info.K * zC + core_info.J * a) / max(norm(core_info.K * zC), 1e-300),
                    B_used = B, B_interface = core_info.B * E1)
        end

        # Love numbers at r_top (R for bc_at = :face): k = Phi(R) - 1 for a unit tidal potential, h = g U(R)
        zR = EN * z[N]      # extrapolate the top node's state to the top radius (face or node)
        kR = Sd[std_of[SΦ]] * zR[loc(SΦ)] - 1
        hR = gtop * Sd[std_of[SU]] * zR[loc(SU)]
        return (k = kR, h = hR, y = y, z = z, faces = faces, slots = sl, types = t, finite = all(isfinite, z[N]),
                r_top = rtop, scale = [Sd[std_of[s]] for s in 1:8], core = core, g_node = gnode, z_C = E1 * z[1],
                z_R = zR)
    end


    # ------------------------------------------------------------------------------------------- verification

    """
        condition_label(c, lo, hi)

    Readable form of a condition, e.g. "U⁻ + 0.2·R⁻ = U⁺", "Y⁻ = 0".
    """
    function condition_label(c, lo::Symbol, hi::Symbol)
        name(s) = ("U", "V", "Φ", "P", "X", "Y", "Ψ", "R")[s]
        side(v, sup) = (ts = [(abs(v[s]) == 1 ? (v[s] < 0 ? "-" : "") : string(round(v[s], sigdigits = 3), "·")) *
                              name(s) * sup for s in 1:8 if v[s] != 0];
                        isempty(ts) ? "0" : join(ts, " + "))
        return side(c[1], "⁻") * " = " * side(c[2], "⁺")
    end

    """
        face_residuals(sol) -> Vector of NamedTuples (i, lo, hi, label, residual, lhs, rhs)

    Residual of every interface condition from the solved node states, with the exact face extrapolations
    y- = C_n z_n and y+ = -D_{n+1} z_{n+1}: |l . y- - g . y+| divided by sum_k (|l_k| max|z_k|- + |g_k| max|z_k|+), the
    maxima taken over the two layers that meet at the face (non-dimensional states). Round-off is ~1e-12 or below.
    """
    function face_residuals(sol)
        N = length(sol.z)
        out = Any[]
        for f in sol.faces
            t = sol.types
            # extent of the layers below and above the face
            lo0 = f.i; while lo0 > 1 && t[lo0-1] === t[f.i]; lo0 -= 1; end
            hi1 = f.i + 1; while hi1 < N && t[hi1+1] === t[f.i+1]; hi1 += 1; end
            slo, shi = sol.slots[f.i], sol.slots[f.i+1]
            zlo = [maximum(abs(sol.z[q][j]) for q in lo0:f.i) for j in eachindex(slo)]
            zhi = [maximum(abs(sol.z[q][j]) for q in f.i+1:hi1) for j in eachindex(shi)]
            ym = f.C * sol.z[f.i]; yp = -(f.D * sol.z[f.i+1])
            for c in vcat(f.up, f.low)
                ℓ, g = c[1][slo], c[2][shi]
                num = abs(sum(ℓ .* ym) - sum(g .* yp))
                den = sum(abs.(ℓ) .* zlo) + sum(abs.(g) .* zhi)
                push!(out, (i = f.i, lo = f.lo, hi = f.hi, label = condition_label(c, f.lo, f.hi),
                            residual = num / max(den, 1e-300), lhs = sum(ℓ .* ym), rhs = sum(g .* yp)))
            end
        end
        return out
    end


    # ------------------------------------------------------------------------------------------------ heating

    """
        column_heating(cs, y, n, ω; solid=:displacement, jump=nothing)

    Dissipation of the column solution `y` (column_solve) per solver cell, from `common.get_heating_profile_nodal`: fields
    `shear`, `bulk`, `darcy`, `drag` (shell-mean densities, W m^-3 per |U_T|^2) and the integrated power `P_total`
    (W per |U_T|^2).

    Strain recovery (`common.strain_parts`): fluid nodes always use the displacement-based strain (y4 / mu_f is ill
    conditioned for a fluid); solid and mush nodes use the same (solid = :displacement) or the stress-based node strain
    (solid = :stress). `jump` (default: every change of node type and every porosity change inside a mush) marks the
    faces across which no displacement difference is taken.
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
    # displacement difference is taken (default: every change of node type and every porosity change inside a mush)
    function heating_masks(cs, solid::Symbol, jump)
        solid in (:displacement, :stress) || error("column_heating: solid must be :displacement or :stress")
        Nc = length(cs.r_c)
        jmp = jump === nothing ? [cs.types[j] !== cs.types[j+1] || (cs.types[j] === :mush && cs.ϕ[j] != cs.ϕ[j+1])
                                  for j in 1:Nc-1] : jump
        return cs.fluid .| (solid == :displacement), cs.types .=== :mush, jmp
    end

    """
        boundary_layer_resolution(cs)

    Poroelastic boundary-layer thickness of every mush node over its cell size (appendix eq:delta):

        delta = sqrt(2 k |M_B| / (eta_l omega)),   1/M_B = S + alpha^2 / beta_d,   beta_d = K_d + 4 mu / 3,
        S = phi / K_l + (alpha - phi) / K_s.

    Darcy displacement R vanishes at an impermeable face and recovers over delta, and at a porosity jump R jumps and
    relaxes over delta on each side. A collapsed mush (K_d = 0, mu ~ G_inf) has |M_B| ~ beta_d and a very thin
    boundary layer (metres for G_inf ~ 1 kPa), which a grid with delta / dr < 1 does not resolve: the energy closure
    then drifts with dr. Returns (`ratio`: per node, NaN for non-mush nodes; `min_ratio`; `n_under`: nodes with
    ratio < 2).
    """
    function boundary_layer_resolution(cs)
        Nc = length(cs.r_c)
        ratio = fill(NaN, Nc)
        for q in 1:Nc
            cs.types[q] === :mush || continue
            βd = cs.κd[q] + 4cs.μ[q] / 3
            S = cs.ϕ[q] / cs.κl[q] + (cs.α[q] - cs.ϕ[q]) / cs.κs[q]
            MB = 1 / (S + cs.α[q]^2 / βd)
            δ = sqrt(2 * cs.k[q] * abs(MB) / (cs.ηl[q] * abs(cs.σ)))
            ratio[q] = δ / (cs.r_grid[q+1] - cs.r_grid[q])
        end
        ok = .!isnan.(ratio)
        return (ratio = ratio, min_ratio = any(ok) ? minimum(ratio[ok]) : Inf, n_under = count(x -> !isnan(x) && x < 2, ratio))
    end

    """
        tidal_work(k, n, r, ω)

    Tidal work per |U_T|^2 [W] for a Love number k of degree n referenced at radius r [m]:
    (2n+1) r w / (8 pi G) (-Im k). Equals the integrated dissipation for an energy-consistent solution.
    """
    tidal_work(k, n, r, ω) = (2n + 1) * r * ω / (8π * G) * (-imag(k))

    """
        solve_column(cs; n=2, solid=:auto, core_route=:auto)

    One frequency of the column model: column_solve, heating and energy closure. solid = :auto evaluates both strain
    recoveries for the solid/mush nodes (fluid nodes are always displacement-based) and keeps the one whose closure
    P_total / tidal work is closer to 1 (both are consistent discretisations of the same dissipation; the closure picks
    the one that resolves the boundary layers of this solution better).

    # Returns
    NamedTuple: `k`, `h`, `sol` (column_solve), `heating` (column_heating of the chosen recovery), `P_tide`,
    `closure`, `closure_disp`, `closure_stress`, `passive` (-sign(ω) Im k >= 0).
    """
    function solve_column(cs; n::Int = 2, solid::Symbol = :auto, core_route::Symbol = :auto)
        sol = column_solve(cs, n; core_route = core_route)
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

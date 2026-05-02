

module solid1d_mush_relax
    
    include("common.jl")
    using .common

    using LinearAlgebra
    import GenericLinearAlgebra
    using DoubleFloats
    using AssociatedLegendrePolynomials    
    using StaticArrays
    using SpecialFunctions
    using SparseArrays
    using Printf

    prec = BigFloat
    precc = Complex{BigFloat}

    const G::prec       = prec(6.6743e-11)       # m^3 kg^-1 s^-2

    clats = 0.0
    lons  = 0.0
    Y     = 0.0
    dYdθ  = 0.0
    dYdϕ  = 0.0
    Z     = 0.0
    X     = 0.0
    res   = 20.0


    """
        resample_profiles(radius, rho, visc, shear, bulk, phi, m_core, dr_min, dr_max)

    Resample the input profiles onto a new grid with `ncalc` points. The new grid is generated using a 
    stretched and refined scheme, which allows for better resolution in regions of interest (e.g., near 
    layer boundaries). 

    # Arguments
    - `radius::Vector{Float64}`           : Original radius profile (layer boundaries).
    - `rho::Vector{Float64}`              : Original density profile (defined at layer centers).
    - `visc::Vector{Float64}`             : Original viscosity profile (defined at layer centers).
    - `shear::Vector{Float64}`            : Original shear modulus profile (defined at layer centers).
    - `bulk::Vector{Float64}`             : Original bulk modulus profile (defined at layer centers).
    - `phi::Vector{Float64}`              : Original phi profile (defined at layer centers).
    - `m_core::Float64`                   : Mass of the core, used for gravity calculations.
    - `Δr_min::Float64`                   : Minimum grid spacing for the new grid.
    - `Δr_max::Float64`                   : Maximum grid spacing for the new grid.

    # Returns
    Tuple of resampled profiles on the new grid:
    - `r_new_b::Vector{prec}`             : New radius profile at layer boundaries.
    - `ρ_new::Vector{prec}`               : New density profile at layer centers.
    - `η_new::Vector{prec}`               : New viscosity profile at layer centers.
    - `μ_new::Vector{prec}`               : New shear modulus profile at layer centers.
    - `κ_new::Vector{prec}`               : New bulk modulus profile at layer centers.
    - `φ_new::Vector{prec}`               : New phi profile at layer centers.
    - `g_new::Vector{prec}`               : New gravity profile at layer centers.
    """ 
    function resample_profiles(radius, rho, visc, shear, bulk_s, bulk_l, bulk_d, alpha, visc_l, phi, k, m_core, dr_min, dr_max)
        # setup grids
        α = log(dr_max / dr_min)

        N = Int(round((radius[end] - radius[1]) / dr_min * α / (exp(α) - 1)))

        # indices i = 1:N
        i = collect(1:N)

        # convert to BigFloat for consistency
        i_bf = BigFloat.(i)
        N_bf = BigFloat(N)

        # compute normalized coordinate (N - i)/(N - 1)
        ξ = (N_bf .- i_bf) ./ (N_bf - 1)

        # compute r_i
        r_new_b = radius[end] .+ (radius[1] - radius[end]) .* (
            (exp.(α .* ξ) .- 1) ./ (exp(α) - 1)
        )

        # cell centers
        r_new_c = 0.5 .* (r_new_b[1:end-1] .+ r_new_b[2:end])

        # obtain new profiles (Constant per original layer)
        ρ_new = similar(rho, N-1)
        η_new = similar(visc, N-1)
        μ_new = similar(shear, N-1)
        κs_new = similar(bulk_s, N-1)
        κl_new = similar(bulk_l, N-1)
        κd_new = similar(bulk_d, N-1)
        α_new = similar(alpha, N-1)
        ηl_new = similar(visc_l, N-1)
        φ_new = similar(phi, N-1)
        k_new = similar(k, N-1)

        for i in 1:N-1
            # find index such that r_b[idx] <= r_new_c[i] < r_b[idx+1]
            idx = searchsortedfirst(radius, r_new_c[i]) - 1
            idx = clamp(idx, 1, length(rho)) # Safety clamp

            ρ_new[i] = rho[idx]
            η_new[i] = visc[idx]
            μ_new[i] = shear[idx]
            κs_new[i] = bulk_s[idx]
            κl_new[i] = bulk_l[idx]
            κd_new[i] = bulk_d[idx]
            α_new[i] = alpha[idx]
            ηl_new[i] = visc_l[idx]
            φ_new[i] = phi[idx]
            k_new[i] = k[idx]
        end

        g_new, M_tot = get_g(r_new_b, ρ_new, m_core) 

        return r_new_b, ρ_new, η_new, μ_new, κs_new, κl_new, κd_new, α_new, ηl_new, φ_new, k_new, g_new, M_tot
    end


    """
        get_g(r, ρ, m_core)

    Compute the radial gravity structure associated with a density profile `r` at intervals given by `r`.

    # Arguments
    - `r::Array{Float64,1}`               : 1D array of layer boundaries. 
    - `ρ::Array{Float64,1}`               : 1D array of layer densities. The length of `ρ` must be equal to the number of columns in `r`.
    - `m_core::Float64`                   : Mass of the core.

    # Returns
    - `g::Array{Float64,1}`               : 1D array of gravity values at the layer boundaries. The dimensions of `g` are the same as `r`.
    - `M_enc::Float64`                    : Total mass enclosed within the outermost layer boundary.
    """
    function get_g(r, ρ, m_core)

        dm = 4.0/3.0 * π .* diff(r.^3) .* ρ

        M_enc = cumsum(dm) .+ m_core
            
        g = G .* M_enc ./ r[2:end].^2

        return g, M_enc[end]
    end


    """
        Ynm(n, m, theta, phi)

    Compute the spherical harmonic Ynm for given n, m, theta, and phi.

    # Arguments
    - `n::Int`                          : Tidal degree.
    - `m::Int`                          : Tidal order.
    - `theta::Array{Float64,1}`         : Array of colatitudes in radians.
    - `phi::Array{Float64,1}`           : Array of longitudes in radians.

    # Returns
    - `Ynm::Array{ComplexF64,2}`        : 2D array of spherical harmonic values for each combination of theta and phi.
    """
    function Ynm(n, m, theta, phi)
        return Plm.(n, m, cos.(theta)) .* exp.(1im * m .* phi)
    end


    """
        define_spherical_grid(res)

    Create the spherical grid of angular resolution `res` in degrees. This is used for 
    numerical integrations over solid angle. A new grid can easily be defined by 
    recalling the function with a new `res`.

    # Arguments
    - `res::Float64`                     : Desired angular resolution in degrees.
    - `n::Int`                           : Tidal degree.
    - `m::Int`                           : Tidal order.

    # Notes
    The grid is internal to solid1d_relax, but can be accessed with 

        solid1d_relax.clats[:] # colatitude grid
        solid1d_relax.lons[:]  # longitude grid
    """
    function define_spherical_grid(res, n, m)
        solid1d_mush_relax.res = res

        # θ and φ grids
        lons = deg2rad.(collect(0:res:360-0.001))'
        clats = deg2rad.(collect(0:res:180))
        clats[1] += 1e-6
        clats[end] -= 1e-6

        # allocate arrays
        solid1d_mush_relax.Y    = zeros(ComplexF64, 1, length(clats), length(lons))
        solid1d_mush_relax.dYdθ = similar(solid1d_mush_relax.Y)
        solid1d_mush_relax.dYdϕ = similar(solid1d_mush_relax.Y)
        solid1d_mush_relax.Z    = similar(solid1d_mush_relax.Y)
        solid1d_mush_relax.X    = similar(solid1d_mush_relax.Y)

        sinθ = sin.(clats)
        cosθ = cos.(clats)
        cotθ = cosθ ./ sinθ
        cscθ = csc.(clats)

        # Normalization factor for spherical harmonics
        norm = sqrt((2*n+1)  * factorial(n-m) / (4π * factorial(n+m)))
        
        i = 1

        # Y
        solid1d_mush_relax.Y[i,:,:] = Ynm(n,m,clats,lons)

        # ∂Y/∂θ
        Pn = Plm.(n, m, cosθ)
        if n > m
            Pn_1 = Plm.(n-1, m, cosθ)
            dPdθ = (n .* cosθ .* Pn .- (n + m) .* Pn_1) ./ (sinθ)
        else
            # m == n -> P_{n-1}^m = 0
            dPdθ = (n .* cosθ .* Pn) ./ (sinθ)
        end
        solid1d_mush_relax.dYdθ[i,:,:] .= dPdθ .* exp.(1im .* m .* lons)

        # ∂Y/∂ϕ
        solid1d_mush_relax.dYdϕ[i,:,:] .= 1im * m .* solid1d_mush_relax.Y[i,:,:]

        # Z = 2 ((1/sinθ) ∂²Y/∂θ∂ϕ - cotθ cscθ ∂Y/∂ϕ)
        solid1d_mush_relax.Z[i,:,:] .= 2 .* (1im * m ./ sinθ .* solid1d_mush_relax.dYdθ[i,:,:] .- cotθ .* cscθ .* solid1d_mush_relax.dYdϕ[i,:,:])

        # X = -2 (cotθ ∂Y/∂θ + csc²θ ∂²Y/∂ϕ²) - n(n+1)) Y
        solid1d_mush_relax.X[i,:,:] .= -2 .* (cotθ .* solid1d_mush_relax.dYdθ[i,:,:] .- cscθ.^2 .* m^2 .* solid1d_mush_relax.Y[i,:,:]) .- n*(n+1) .* solid1d_mush_relax.Y[i,:,:]

        # Normalize
        solid1d_mush_relax.Y[i,:,:]    .*= norm
        solid1d_mush_relax.dYdθ[i,:,:] .*= norm
        solid1d_mush_relax.dYdϕ[i,:,:] .*= norm
        solid1d_mush_relax.Z[i,:,:]    .*= norm
        solid1d_mush_relax.X[i,:,:]    .*= norm

        # save grids
        solid1d_mush_relax.clats = clats
        solid1d_mush_relax.lons  = lons
    end

    
    """
        solve_radial_system(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n, ρ_core; core="liquid")

    Solve the radial system of ODEs for the solid-body problem using a relaxation method. This function 
    implements the forward-backward relaxation scheme described in the main text of N. Kobayashi (2006).
    
    # Arguments
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `ρₗ::Vector{prec}`                : Vector of liquid densities at the layer centers.
    - `Kl::Vector{prec}`                : Vector of liquid bulk moduli at the layer centers.
    - `Kd::Vector{prec}`                : Vector of drained bulk moduli at the layer centers.
    - `α::Vector{prec}`                 : Vector of Biot coefficients at the layer centers.
    - `ηₗ::Vector{prec}`               : Vector of liquid viscosities at the layer centers.
    - `ϕ::Vector{prec}`                 : Vector of porosities at the layer centers.
    - `k::Vector{prec}`                 : Vector of permeabilities at the layer centers.
    - `n::Int`                          : Tidal degree.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `M_tot::prec`                     : Total mass of the body, used for non-dimensionalization.

    # Keyword Arguments
    - `core::String="liquid"            : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.

    # Returns
    - `y_t::Vector{precc}`              : Vector of length 8 representing the tidal solution at the top of the mantle. This includes the displacements, stresses, and potential at the surface.
    - `y_l::Vector{precc}`              : Vector of length 8 representing the load solution at the top of the mantle. This includes the displacements, stresses, and potential at the surface.
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `S::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the normalization.
    - `transitions::Vector{Int}`        : Indices of the interface layers (the ones closer to the core and surface).
    """
    function solve_radial_system(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n,
                                ρ_core, M_tot; core="liquid")

        # Define the ordering of the solution vector components for the 6x6 and 8x8 cases
        Y6 = [1,2,4,5,3,6]
        Y8 = [1,2,5,6,3,7,4,8]

        # 1. Identify Regions and Transitions
        # We define a boolean mask where true = mushy, false = solid
        is_mush = k .> 0
        
        # Detect indices where the state changes
        transitions = findall(diff(is_mush) .!= 0)
        
        # 2. Duplicate nodes at transition points for the relaxation scheme
        # We iterate backwards to keep indices valid during insertion
        new_r, new_ρ, new_g, new_μ, new_K, new_ρₗ, new_Kl, new_Kd, new_α, new_ηₗ, new_ϕ, new_k = 
            copy(r), copy(ρ), copy(g), copy(μ), copy(K), copy(ρₗ), copy(Kl), copy(Kd), copy(α), copy(ηₗ), copy(ϕ), copy(k)

        for trans_idx in reverse(transitions)
            for v in (new_r, new_ρ, new_g, new_μ, new_K, new_ρₗ, new_Kl, new_Kd, new_α, new_ηₗ, new_ϕ, new_k)
                insert!(v, trans_idx + 1, v[trans_idx])
            end
        end

        Nr = length(new_r)

        new_is_mush = new_k .> 0

        # 3. Dynamic Scaling (needs to be checked for consistency with the non-dimensionalization in get_scales)
        R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(1., 1., 1.; Y=Y8)
        ωs = ω * s0
        rs, ρs, gs, μs, Ks = new_r./R0, new_ρ./ρ0, new_g./g0, new_μ./μ0, new_K./μ0
        ρₗs, Kls, Kds, ηₗs, ks = new_ρₗ./ρ0, new_Kl./μ0, new_Kd./μ0, new_ηₗ./(μ0*s0), new_k./R0^2

        # 4. Initialize Matrices
        R = [Matrix{precc}(I, 8, 8) for _ in 1:Nr]
        B = [zeros(precc, 8, 1) for _ in 1:Nr]
        idx_6 = [1, 2, 3, 5, 6, 7]
        R6_view = [view(R[i], idx_6, idx_6) for i in 1:Nr]
        R8_view = [view(R[i], 1:8, 1:8) for i in 1:Nr]

        # 5. Adaptive Propagation Loop
        @debug("\n--- Adaptive Solver Plan ---")
        @debug("Total Nodes (Nr): $Nr")
        @debug("Transitions at indices: $transitions")

        curr_idx = 1
        # Step 1: Core
        if !new_is_mush[1]
            @debug("STEP 1: [Core Boundary] Solid | Indices: (1, 2)")
            C1l, D2l = core_boundary(R6_view, (1, 2), rs, ρs, gs, μs, Ks, ωs, ρ_core/ρ0, core, n; G0=G0, Y=Y6)
            curr_idx = 2
        else
            @debug("STEP 1: [Core Boundary] Mushy | Indices: (1, 2)")
            C1l, D2l = core_boundary_mush(R8_view, (1, 2), rs, ρs, gs, μs, Ks, ωs, ρₗs, Kls, Kds, new_α, ηₗs, new_ϕ, ks, ρ_core/ρ0, core, n; G0=G0, Y=Y8)
            curr_idx = 2
        end

        # Step 2: Propagation and Jumps
        while curr_idx < Nr
            next_change = findnext(x -> x != new_is_mush[curr_idx], new_is_mush, curr_idx)
            segment_end = (next_change === nothing) ? Nr : next_change - 1
            
            # Safety check for empty ranges
            if segment_end > curr_idx
                if !new_is_mush[curr_idx]
                    @debug("STEP: [Propagate Solid] | Range: ($curr_idx, $segment_end)")
                    C1l, D2l = propagate_solid(R6_view, C1l, D2l, (curr_idx, segment_end-1), 
                                            rs, ρs, gs, μs, Ks, ωs, n; G0=G0, Y=Y6)
                else
                    @debug("STEP: [Propagate Mushy] | Range: ($curr_idx, $segment_end)")
                    C1l, D2l = propagate_mush(R8_view, C1l, D2l, (curr_idx, segment_end-1), 
                                            rs, ρs, gs, μs, Ks, ωs, ρₗs, Kls, Kds, new_α, ηₗs, new_ϕ, ks, n; G0=G0, Y=Y8)
                end
            end

            if next_change !== nothing
                # The jump occurs at the duplicated node
                trans_range = (next_change - 1, next_change + 1)
                if !new_is_mush[curr_idx] && new_is_mush[next_change]
                    @debug("STEP: [Interface Jump] Solid -> Mushy | Range: $trans_range")
                    C1l, D2l = interface_solid_mush(R8_view, C1l, D2l, trans_range; Y=Y8)
                else
                    @debug("STEP: [Interface Jump] Mushy -> Solid | Range: $trans_range")
                    C1l, D2l = interface_mush_solid(R8_view, C1l, D2l, trans_range; Y=Y8)
                end
                curr_idx = next_change + 1
            else
                curr_idx = Nr
            end
        end

        # Step 3: Surface
        if !new_is_mush[Nr]
            @debug("STEP: [Surface Boundary] Solid | Indices: ($(Nr-1), $Nr)")
            y_t, y_l = surface_boundary(R6_view, C1l, D2l, (Nr-1, Nr), rs, ρs, gs, μs, Ks, ωs, n; G0=G0, Y=Y6)
        else
            @debug("STEP: [Surface Boundary] Mushy | Indices: ($(Nr-1), $Nr)")
            y_t, y_l = surface_boundary_mush(R8_view, C1l, D2l, (Nr-1, Nr), rs, ρs, gs, μs, Ks, ωs, n; G0=G0, Y=Y8)
        end
        @debug("--- Solver Complete ---\n")

        return y_t, y_l, R, S, transitions 
    end


    """
        interface_mush_solid(R8, Cn_l, Dnp_l, ids; Y=[1,2,3,4,5,6,7,8])

    Perform the forward-backward relaxation step at the interface between the mushy layer and the solid layer. This 
    function implements the recursion described in N. Kobayashi (2007) for the transition from the 8x8 system to the 6x6 system.

    # Arguments
    - `R8::Vector{Matrix{precc}}`       : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `Cn_l::Matrix{precc}`             : 3x6 matrix representing the "stored" lower half of the Cn matrix from the previous step.
    - `Dnp_l::Matrix{precc}`            : 3x6 matrix representing the "stored" lower half of the Dnp matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.

    # keyword arguments
    - `Y::Vector{Int}=1:8`              : Vector of column indices corresponding to the 6x6 system variables in the 8x8 system. This allows for

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 3x6 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 3x6 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.
    """
    function interface_mush_solid(R8, Cn_l, Dnp_l, ids; Y=[1,2,3,4,5,6,7,8])

        start_id, end_id = ids

        I8  = Matrix{precc}(I, 8, 8)

        Cn  = I8
        Dnp = -I8
        Dnp[4, 4] = 0.0 
        Dnp[8, 8] = 0.0

        target_cols = [Y[1], Y[2], Y[3], Y[4], Y[5], Y[6]]
        # 1. Use the "stored" lower halves from the previous step 
        # to fill the upper blocks of P and S.
        Pn_u = Cn_l
        Sn_u = Dnp_l
        Qn_u = zeros(precc, 4, 8)

        # 2. Get the upper halves of the NEWLY calculated Cn and Dnp
        Cn_u  = Cn[1:4, :]
        Dnp_u = Dnp[1:4, :]

        # 3. Build the 8x8 blocks
        Pn = [Pn_u; zeros(precc, 4, 8)]
        Sn = [Sn_u; Cn_u]
        Qn = [Qn_u; Dnp_u]

        Kn = zeros(precc, 8, 8)
        Kn[8, 8] = 1.0

        # 4. Perform recursion
        Xn = Pn * R8[start_id] + Sn + Kn

        R_ifc = - pinv(Xn) * Qn

        R8[start_id+1] .= R_ifc 

        # 5. Update the "stored" lower halves for the next iteration
        Cn_l  = Cn[5:7, target_cols]
        Dnp_l = Dnp[5:7, target_cols]

        return Cn_l, Dnp_l
    end


    """
        interface_solid_mush(R, Cn_l, Dnp_l, ids; Y=[1,2,3,4,5,6,7,8])

    Perform the forward-backward relaxation step at the interface between the solid layer and the mushy layer. This 
    function implements the recursion described in N. Kobayashi (2007) for the transition from the 6x6 system to the 8x8 system.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `Cn_l::Matrix{precc}`             : 3x6 matrix representing the "stored" lower half of the Cn matrix from the previous step.
    - `Dnp_l::Matrix{precc}`            : 3x6 matrix representing the "stored" lower half of the Dnp matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.

    # keyword arguments
    - `Y::Vector{Int}=1:8`              : Vector of column indices corresponding to the 6x6 system variables in the 8x8 system. This allows for

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 3x6 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 3x6 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.    
    """
    function interface_solid_mush(R, Cn_l, Dnp_l, ids; Y=[1,2,3,4,5,6,7,8])

        start_id, end_id = ids

        # target_cols
        target_cols = [Y[1], Y[2], Y[3], Y[4], Y[5], Y[6]]
        
        # expand the incoming 3x6 Solid Lower blocks to 4x8
        Pn_u = zeros(precc, 4, 8)
        Pn_u[1:3, target_cols] .= Cn_l
        
        Sn_u = zeros(precc, 4, 8)
        Sn_u[1:3, target_cols] .= Dnp_l
                
        # current Layer (Porous side)
        # treat the interface as an infinitesimal jump where C = I, D = -I
        I8 = Matrix{precc}(I, 8, 8)
        Cn_curr = I8
        Cn_curr[4, 4] = 0.0 
        Cn_curr[8, 8] = 0.0
        Dnp_curr = -I8
        
        Cn_u = Cn_curr[1:4, :]
        Dnp_u = Dnp_curr[1:4, :]
        
        # assemble 8x8 system
        Pn = [Pn_u; zeros(precc, 4, 8)]
        Sn = [Sn_u; Cn_u]
        Qn = [zeros(precc, 4, 8); Dnp_u]

        # introduce some pore pressure at the boundary to drive the solution in the mushy layer
        Kn = zeros(precc, 8, 8)
        Kn[8, 8] = 1.0

        # solve the jump
        Xn = Pn * R[start_id] + Sn + Kn
        R_ifc = - pinv(Xn) * Qn

        R[start_id+1] .= R_ifc
        
        # pass the Lower halves of the Porous Identity to the next propagator
        Cn_l_new = Cn_curr[5:8, :]
        Dnp_l_new = Dnp_curr[5:8, :]
        
        return Cn_l_new, Dnp_l_new
    end


    """
        core_boundary(R, ids, r, ρ, g, μ, K, ω, ρ_core, core, n; G0=1, Y=[1,2,3,4,5,6])

    Perform the forward-backward relaxation step at the core boundary. This function implements the recursion described 
    in N. Kobayashi (2007) for the initial step of the relaxation scheme, where we apply the core boundary condition and 
    get the first solution for the first layer above the core.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `core::String`                    : Type of core boundary condition to apply.
    - `n::Int`                          : Tidal degree.

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.
    - `Y::Vector{Int}=[1,2,3,4,5,6]`   : Ordering of the solution vector components.

    # Returns
    - `C1l::Matrix{precc}`              : 3x6 matrix representing the "stored" lower half of the C1 matrix for the next iteration.
    - `D2l::Matrix{precc}`              : 3x6 matrix representing the "stored" lower half of the D2 matrix for the next iteration.
    """
    function core_boundary(R, ids, r, ρ, g, μ, K, ω, ρ_core, core, n; G0=1, Y=[1,2,3,4,5,6])

        start_id, end_id = ids

        # boundary conditions
        B1 = get_core_bc!(ω, r[start_id], ρ_core, g[start_id], μ[start_id], K[start_id], core, n; G0=G0, M=6, N=3)        
        
        # first layer (n = 1)
        dr = r[end_id] - r[start_id]

        A1 = get_A(ω, r[start_id], ρ[start_id], g[start_id], μ[start_id], K[start_id], n; G0=G0, Y=Y)
        A2 = get_A(ω, r[end_id], ρ[end_id], g[end_id], μ[end_id], K[end_id], n; G0=G0, Y=Y)

        I6 = Matrix{precc}(I, 6, 6)

        C1 =  I6 + 0.5 * dr * A1
        D2 = -I6 + 0.5 * dr * A2

        # split matrices
        C1u, C1l = C1[1:3, :], C1[4:6, :]
        D2u, D2l = D2[1:3, :], D2[4:6, :]

        # build S1 and Q1
        S1 = [B1; C1u]              # 6×6
        Q1 = [zeros(3,6); D2u]      # 6×6

        # initial recursion
        R[start_id] .= -S1 \ Q1

        return C1l, D2l
    end


    """
        core_boundary_mush(R, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, ρ_core, core, n; G0=1, Y=[1,2,3,4,5,6,7,8])

    Perform the forward-backward relaxation step at the core boundary for the two-phase problem. This function implements 
    the recursion described in N. Kobayashi (2007) for the initial step of the relaxation scheme, where we apply the core 
    boundary condition and get the first solution for the first layer above the core, but now using the full 8x8 system 
    that includes the porous layer components.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `ρₗ::Vector{prec}`                : Vector of liquid densities at the layer centers.
    - `Kl::Vector{prec}`                : Vector of liquid bulk moduli at the layer centers.
    - `Kd::Vector{prec}`                : Vector of drained bulk moduli at the layer centers.
    - `α::Vector{prec}`                 : Vector of Biot coefficients at the layer centers.
    - `ηₗ::Vector{prec}`               : Vector of liquid viscosities at the layer centers.
    - `ϕ::Vector{prec}`                 : Vector of porosities at the layer centers.
    - `k::Vector{prec}`                 : Vector of permeabilities at the layer centers.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `core::String`                    : Type of core boundary condition to apply.
    - `n::Int`                          : Tidal degree.

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.

    # Returns
    - `C1l::Matrix{precc}`              : 4x8 matrix representing the "stored" lower half of the C1 matrix for the next iteration.
    - `D2l::Matrix{precc}`             : 4x8 matrix representing the "stored" lower half of the D2 matrix for the next iteration.
    """
    function core_boundary_mush(R, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, ρ_core, core, n; G0=1, Y=[1,2,3,4,5,6,7,8])

        start_id, end_id = ids

        # boundary conditions
        B1 = get_core_bc!(ω, r[start_id], ρ_core, g[start_id], μ[start_id], K[start_id], core, n; G0=G0, Y=Y)     
        
        # first layer (n = 1)
        dr = r[end_id] - r[start_id]

        A1 = get_A(ω, r[start_id], ρ[start_id], g[start_id], μ[start_id], K[start_id],
                    ρₗ[start_id], Kl[start_id], Kd[start_id], α[start_id], ηₗ[start_id], ϕ[start_id], k[start_id], n; G0=G0, Y=Y)

        A2 = get_A(ω, r[end_id], ρ[end_id], g[end_id], μ[end_id], K[end_id],
                ρₗ[end_id], Kl[end_id], Kd[end_id], α[end_id], ηₗ[end_id], ϕ[end_id], k[end_id], n; G0=G0, Y=Y)

        I8 = Matrix{precc}(I, 8, 8)

        C1 =  I8 + 0.5 * dr * A1
        D2 = -I8 + 0.5 * dr * A2

        # split matrices
        C1u, C1l = C1[1:4, :], C1[5:8, :]
        D2u, D2l = D2[1:4, :], D2[5:8, :]

        # build S1 and Q1
        S1 = [B1; C1u]              # 8×8
        Q1 = [zeros(4,8); D2u]      # 8×8

        # initial recursion
        R[start_id] .= -pinv(S1) * Q1

        return C1l, D2l
    end

    
    """
        propagate_solid(R, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, n; G0=1, Y=[1,2,3,4,5,6])

    Perform the forward-backward relaxation step for the solid propagation segments. This function implements the 
    recursion described in N. Kobayashi (2007) for the segments of the radial grid that correspond to the solid 
    layers, where we propagate the solution up to the surface using the 6x6 system of equations.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `Cn_l::Matrix{precc}`             : 3x6 matrix representing the "stored" lower half of the Cn matrix from the previous step.
    - `Dnp_l::Matrix{precc}`            : 3x6 matrix representing the "stored" lower half of the Dnp matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `n::Int`                          : Tidal degree.

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.
    - `Y::Vector{Int}=[1,2,3,4,5,6]`    : Ordering of the solution vector components.

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 3x6 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 3x6 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.
    """
    function propagate_solid(R, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, n; G0=1, Y=[1,2,3,4,5,6])

        start_id, end_id = ids

        I6 = Matrix{precc}(I, 6, 6)

        Cn_u = zeros(3,6)
        Dnp_u = zeros(3,6)

        # forward recursion
        for i in start_id:end_id

            dr = r[i+1] - r[i]

            # Calculate A at current and next step
            A_n  = get_A(ω, r[i],   ρ[i],   g[i],   μ[i],   K[i],   n; G0=G0, Y=Y)
            A_np = get_A(ω, r[i+1], ρ[i+1], g[i+1], μ[i+1], K[i+1], n; G0=G0, Y=Y)

            Cn  =  I6 + 0.5 * dr * A_n
            Dnp = -I6 + 0.5 * dr * A_np

            # 1. Use the "stored" lower halves from the previous step 
            # to fill the upper blocks of P and S.
            Pn_u = Cn_l
            Sn_u = Dnp_l
            Qn_u = zeros(precc, 3, 6)

            # 2. Get the upper halves of the NEWLY calculated Cn and Dnp
            Cn_u  = Cn[1:3, :]
            Dnp_u = Dnp[1:3, :]

            # 3. Build the 6x6 blocks
            Pn = [Pn_u; zeros(precc, 3, 6)]
            Sn = [Sn_u; Cn_u]
            Qn = [Qn_u; Dnp_u]

            # 4. Perform recursion
            Xn = Pn * R[i-1] + Sn
            
            if i == start_id
                # For the first step into the mush, we may need to use pinv if the system is not yet fully constrained by the solid boundary conditions.
                R[i] .= - pinv(Xn) * Qn
            else
                R[i] .= -Xn \ Qn
            end
            
            # 5. Update the "stored" lower halves for the next iteration
            Cn_l  = Cn[4:6, :]
            Dnp_l = Dnp[4:6, :]
        end

        return Cn_l, Dnp_l
    end

    
    """
        propagate_mush(R, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1, Y=[1,2,3,4,5,6,7,8])

    Perform the forward-backward relaxation step for the mushy layer propagation segment. This function implements the
    recursion described in N. Kobayashi (2007) for the segment of the radial grid that corresponds to the mushy layer, 
    where we propagate the solution up to the surface using the full 8x8 system of equations that includes the porous 
    layer components.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `Cn_l::Matrix{precc}`             : 4x8 matrix representing the "stored" lower half of the Cn matrix from the previous step.
    - `Dnp_l::Matrix{precc}`            : 4x8 matrix representing the "stored" lower half of the Dnp matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `ρₗ::Vector{prec}`                : Vector of liquid densities at the layer centers.
    - `Kl::Vector{prec}`                : Vector of liquid bulk moduli at the layer centers.
    - `Kd::Vector{prec}`                : Vector of drained bulk moduli at the layer centers.
    - `α::Vector{prec}`                 : Vector of Biot coefficients at the layer centers.
    - `ηₗ::Vector{prec}`               : Vector of liquid viscosities at the layer centers.
    - `ϕ::Vector{prec}`                 : Vector of porosities at the layer centers.
    - `k::Vector{prec}`                 : Vector of permeabilities at the layer centers.
    - `n::Int`                          : Tidal degree. 

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.  
    - `Y::Vector{Int}=[1,2,3,4,5,6,7,8]` : Indices for the state variables (default is for standard case).

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 4x8 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 4x8 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.
    """
    function propagate_mush(R, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1, Y=[1,2,3,4,5,6,7,8])

        start_id, end_id = ids

        I8 = Matrix{precc}(I, 8, 8)

        Cn_u = zeros(4,8)
        Dnp_u = zeros(4,8)

        # forward recursion
        for i in start_id:end_id

            dr = r[i+1] - r[i]

            # Calculate A at current and next step
            A_n  = get_A(ω, r[i],   ρ[i],   g[i],   μ[i],   K[i],
                            ρₗ[i], Kl[i], Kd[i], α[i], ηₗ[i], ϕ[i], k[i], n; G0=G0, Y=Y)

            A_np = get_A(ω, r[i+1], ρ[i+1], g[i+1], μ[i+1], K[i+1],
                        ρₗ[i+1], Kl[i+1], Kd[i+1], α[i+1], ηₗ[i+1], ϕ[i+1], k[i+1], n; G0=G0, Y=Y)

            Cn  =  I8 + 0.5 * dr * A_n
            Dnp = -I8 + 0.5 * dr * A_np

            # 1. Use the "stored" lower halves from the previous step 
            # to fill the upper blocks of P and S.
            Pn_u = Cn_l
            Sn_u = Dnp_l
            Qn_u = zeros(precc, 4, 8)

            # 2. Get the upper halves of the NEWLY calculated Cn and Dnp
            Cn_u  = Cn[1:4, :]
            Dnp_u = Dnp[1:4, :]

            # 3. Build the 6x6 blocks
            Pn = [Pn_u; zeros(precc, 4, 8)]
            Sn = [Sn_u; Cn_u]
            Qn = [Qn_u; Dnp_u]

            # 4. Perform recursion
            Xn = Pn * R[i-1] + Sn

            if i == start_id
                # For the first step into the mush, we may need to use pinv if the system is not yet fully constrained by the solid boundary conditions.
                R[i] .= -pinv(Xn) * Qn
            else
                R[i] .= -Xn \ Qn
            end

            # 5. Update the "stored" lower halves for the next iteration
            Cn_l  = Cn[5:8, :]
            Dnp_l = Dnp[5:8, :]
        end

        return Cn_l, Dnp_l
    end


    """
        surface_boundary(R, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n)

    Perform the forward-backward relaxation step at the surface boundary. This function implements the recursion described 
    in N. Kobayashi (2007) for the final step of the relaxation scheme, where we apply the surface boundary condition and 
    solve for the final solution at the surface, using the 6x6 system of equations.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `CNm_l::Matrix{precc}`            : 3x6 matrix representing the "stored" lower half of the CNm matrix from the previous step.
    - `DN_l::Matrix{precc}`             : 3x6 matrix representing the "stored" lower half of the DN matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `n::Int`                          : Tidal degree.

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.

    # Returns
    - `y_t::Matrix{precc}`              : 6x1 matrix representing the solution at the surface for the tidal problem.
    - `y_l::Matrix{precc}`              : 6x1 matrix representing the solution at the surface for the load problem.
    """
    function surface_boundary(R, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n; G0=1, Y=Y)

        start_id, end_id = ids

        # tidal surface boundary condition
        BN_t, b_t = get_surface_bc!(r[end], g[end], n, 1, 0, 0, 0; G0=G0, Y=Y)
        # load surface boundary condition
        BN_l, b_l = get_surface_bc!(r[end], g[end], n, 0, 1, 0, 0; G0=G0, Y=Y)

        PN = [CNm_l; zeros(3,6)]
        SN_t = [DN_l; BN_t]
        SN_l = [DN_l; BN_l]

        XN_t = PN * R[start_id] + SN_t
        XN_l = PN * R[start_id] + SN_l

        # solve outer  (tides)
        y_t = XN_t \ b_t
        # solve outer (load)
        y_l = XN_l \ b_l

        return y_t, y_l

    end


    """
        surface_boundary_mush(R, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n)

    Perform the forward-backward relaxation step at the surface boundary for the two-phase problem. This function implements 
    the recursion described in N. Kobayashi (2007) for the final step of the relaxation scheme, where we apply the surface 
    boundary condition and solve for the final solution at the surface, using the full 8x8 system of equations that includes 
    the porous layer components.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `CNm_l::Matrix{precc}`            : 4x8 matrix representing the "stored" lower half of the CNm matrix from the previous step.
    - `DN_l::Matrix{precc}`             : 4x8 matrix representing the "stored" lower half of the DN matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `n::Int`                          : Tidal degree.

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.

    # Returns
    - `y_t::Matrix{precc}`              : 8x1 matrix representing the solution at the surface for the tidal problem.
    - `y_l::Matrix{precc}`              : 8x1 matrix representing the solution at the surface for the load problem.
    """
    function surface_boundary_mush(R, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n; G0=1, Y=[1,2,3,4,5,6,7,8])

        start_id, end_id = ids

        # tidal surface boundary condition
        BN_t, b_t = get_surface_bc!(r[end], g[end], n, 1, 0, 0, 0; Y=Y)
        # load surface boundary condition
        BN_l, b_l = get_surface_bc!(r[end], g[end], n, 0, 1, 0, 0; Y=Y)

        PN = [CNm_l; zeros(4,8)]
        SN_t = [DN_l; BN_t]
        SN_l = [DN_l; BN_l]

        XN_t = PN * R[start_id] + SN_t
        XN_l = PN * R[start_id] + SN_l

        # solve outer  (tides)
        y_t = pinv(XN_t) * b_t
        # solve outer (load)
        y_l = pinv(XN_l) * b_l

        return y_t, y_l

    end


    """
        get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1, Y=[1,2,3,4,5,6])

    Get the core boundary condition matrix `B` for the solid-body problem. The core boundary 
    conditions are derived from the requirement that the radial stress at the core-mantle 
    boundary must balance the tidal potential, and that the tangential stresses must vanish.

    # Arguments
    - `ω::Float64`                       : Forcing frequency.
    - `r::prec`                          : Radial position of the core-mantle boundary.
    - `ρ::prec`                          : Average core density.
    - `g::prec`                          : Gravity at the core-mantle boundary.
    - `μ::precc`                         : Average core complex shear modulus.
    - `K::prec`                          : Average core bulk modulus.
    - `type::String`                     : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.
    - `Y::Vector{Int}=[1,2,3,4,5,6]`    : Indices for the state variables (default is for standard case).

    # Returns
    - `B::Array{precc,2}`                : 3x6 matrix representing the linear relationship between the state variables at the core and the boundary conditions.
    """
    function get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1, Y=[1,2,3,4,5,6])
        
        M = length(Y)
        N = Int(M / 2)

        # 1. Get the Initial Conditions matrix
        Ic = get_Ic(ω, r, ρ, g, μ, K, type, n; G0=G0, Y=Y)

        # 2. Define indices based on dimensionality
        # If M=8 (Mushing/Hay 2025):  U=1, V=2, phi=5, P=7 | X=3, Y=4, psi=6, R=8
        # If M=6 (Standard/Takeuchi): U=1, V=2, phi=5      | X=3, Y=4, psi=6
        if M == 8
            idx_u = [Y[1], Y[2], Y[5], Y[7]]
            idx_s = [Y[3], Y[4], Y[6], Y[8]]
        else
            idx_u = [Y[1], Y[2], Y[5]]
            idx_s = [Y[3], Y[4], Y[6]]
        end

        # 3. Partition and calculate the boundary condition coefficients
        Mu = Ic[idx_u, :]
        Ms = Ic[idx_s, :]
        
        # Equation 91: b = -Mu * Ms⁻¹
        b = -Mu * pinv(Ms)

        # 4. Construct the NxM B matrix
        T = eltype(b)
        B = zeros(T, N, M)

        for i in 1:N
            B[i, idx_u[i]] = 1.0
            for j in 1:length(idx_s)
                B[i, idx_s[j]] = b[i, j]
            end
        end

        return B
    end


    """
        get_surface_bc!(R, g, n, U, U_prime, tau, P; G0=1, Y=[1,2,3,4,5,6,7,8])

    Get the surface boundary condition vector `b` and matrix `BN` for the solid-body problem. The surface 
    boundary conditions are determined by setting, respectively (U, U', tau, P) to (1,0,0,0) for tidal Love 
    number and (0,1,0,0) for load Love number in system.

    https://hal.science/hal-03421553/document

    # Arguments
    - `R::prec`                          : Planetary radius, used for surface boundary conditions.
    - `g::prec`                          : Gravity at the surface, used for surface boundary conditions.
    - `n::Int`                           : Tidal degree.
    - `U::Int`                           : Tidal potential at the surface.
    - `U_prime::Int`                     : Radial derivative of the tidal potential at the surface.
    - `tau::Int`                         : Tangential tidal stress at the surface.
    - `P::Int`                           : Surface mass load at the surface.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.
    - `Y::Vector{Int}=[1,2,3,4,5,6,7,8]` : Indices for the state variables (default is for standard case).

    # Returns
    - `B::Array{precc,2}`                : 3x6 matrix representing the linear relationship between the state variables at the surface and the boundary conditions.
    - `b::Vector{precc}`                 : Vector of length 6 representing the inhomogeneous part of the surface boundary conditions.
    """
    function get_surface_bc!(R, g, n, U, U_prime, tau, P; G0=1, Y=[1,2,3,4,5,6,7,8])
        
        M = length(Y)
        N = Int(M / 2)

        # Define surface mass load (zeta) based on Farrell/Longman relation
        zeta = ((2 * n + 1) / (4 * pi * G/G0 * R)) * U_prime

        # b vector (Right Hand Side of the B*y = b system)
        b = zeros(precc, M) 
        
        if M == 8
            # radial Stress y3
            b[Y[3]] = -g * zeta * (G/G0) / R - P
            
            # tangential Stress y4
            b[Y[4]] = tau
            
            # gravitational potential boundary
            b[Y[6]] = ((2 * n + 1) / R) * U + 4 * pi * G/G0 * zeta

            # darcy flux boundary
            b[Y[8]] = 0
        elseif M == 6
            # radial Stress y3
            b[Y[3]] = -g * zeta * (G/G0) / R - P
            
            # tangential Stress y4
            b[Y[4]] = tau
            
            # gravitational potential boundary
            b[Y[6]] = ((2 * n + 1) / R) * U - 4 * pi * G/G0 * zeta
        else
            error("Unsupported M value. M should be either 6 or 8.")
        end
        
        # construct the 4x8 B matrix
        # this matrix extracts y3, y4, and the combination for y6
        B = zeros(precc, N, M)

        if M == 8
            # stress components
            B[1, Y[3]] = 1.0  # radial stress y3
            B[2, Y[4]] = 1.0  # tangential stress y4
            
            # potential component
            B[3, Y[5]] = (n + 1) / R
            B[3, Y[6]] = 1.0
            B[4, Y[8]] = 1.0
        elseif M == 6
            # stress components
            B[1, Y[3]] = 1.0  # radial stress y3
            B[2, Y[4]] = 1.0  # tangential stress y4
            # potential component
            B[3, Y[5]] = (n + 1) / R
            B[3, Y[6]] = 1.0        
        end

        return B, b
    end


    """
        compute_y(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n, ρ_core; core="liquid")

    Compute the solution `y` to the solid-body problem using a relaxation method. This function performs the 
    forward-backward relaxation scheme described in the main text of N. Kobayashi (2006), where we first solve 
    the radial system of ODEs to obtain the solution at the surface, and then perform back-substitution to 
    compute the solution at all interior points. 

    # Arguments
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `ρₗ::Vector{prec}`                : Vector of liquid densities at the layer centers.
    - `Kl::Vector{prec}`                : Vector of liquid bulk moduli at the layer centers.
    - `Kd::Vector{prec}`                : Vector of drained bulk moduli at the layer centers.
    - `α::Vector{prec}`                 : Vector of Biot coefficients at the layer centers.
    - `ηₗ::Vector{prec}`               : Vector of liquid viscosities at the layer centers.
    - `ϕ::Vector{prec}`                 : Vector of porosities at the layer centers.
    - `k::Vector{prec}`                 : Vector of permeabilities at the layer centers.
    - `n::Int`                          : Tidal degree.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `M_tot::prec`                     : Total mass of the planet.

    # Keyword Arguments
    - `core::String="liquid"            : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.

    # Returns
    - `y::Matrix{precc}`                : 6xN matrix of the solution at all radial grid points, where N is the number of radial layers. Each column corresponds to a radial grid point, and each row corresponds to a state variable (displacements, stresses, potential).
    """    
    function compute_y(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n, ρ_core, M_tot; core="liquid")

        # solve radial system to get surface solution and recursion matrices
        yN_t, yN_l, R, S, transitions = solve_radial_system(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n, ρ_core, M_tot; core=core)

        Nr = length(r) + length(transitions) 
        T = eltype(yN_t)
        M = length(yN_t)

        # allocate as 8 x N matrix 
        y_t = zeros(T, 8, Nr)
        y_l = zeros(T, 8, 1)

        # solve outer  (tides)
        y_t[1:M, Nr] = yN_t
        # solve outer (load)
        y_l[1:M, 1] = yN_l

        # back-substitution
        for i in Nr-1:-1:1
            # if this is not a transition point, we perform the recursion step as normal
            y_t[:, i] = R[i] * y_t[:, i+1]        
        end

        # keep only indices NOT in the transition list
        mask = filter(i -> !(i in transitions), 1:size(y_t, 2))

        # Apply the mask to the columns
        y_t = y_t[:, mask]

        # apply scaling to get dimensional solution
        # for i in 1:Nr-2
        for i in 1:Nr-length(transitions)
            y_t[:, i] = S * y_t[:, i] 
        end

        y_l .= S * y_l 

        return y_t, y_l
    end


    """
        compute_strain_ten!(ϵ, y, n, rr, ρr, gr, μr, Ksr, ω, ρlr, Klr, Kdr, αr, ηlr, ϕr, kr)

    Calculate the strain tensor ϵ at a particular radial level. 

    # Arguments
    - `ϵ::Array{ComplexF64,4}`          : 4D array to store the strain tensor at a particular radial level, with dimensions corresponding to latitude, longitude, and the 6 independent components of the strain tensor.
    - `y::Array{precc,1}`               : 1D array of the solution vector y at a particular radial level, with 6 components.
    - `n::Int`                          : Tidal degree.
    - `rr::prec`                        : Radius at which to compute the strain tensor.
    - `ρr::prec`                        : Density at radius rr.
    - `gr::prec`                        : Gravity at radius rr.
    - `μr::prec`                        : Shear modulus at radius rr.
    - `Ksr::prec`                       : Bulk modulus at radius rr.
    - `ω::prec`                         : Forcing frequency.
    - `ρlr::prec`                       : Liquid density at radius rr.
    - `Klr::prec`                       : Liquid bulk modulus at radius rr.
    - `Kdr::prec`                       : Drained bulk modulus at radius rr.
    - `αr::prec`                        : Biot coefficient at radius rr.
    - `ηlr::prec`                       : Liquid viscosity at radius rr.
    - `ϕr::prec`                        : Porosity at radius rr.
    - `kr::prec`                        : Permeability at radius rr.
    """
    function compute_strain_ten!(ϵ, y, n, rr, ρr, gr, μr, Ksr, ω, ρlr, Klr, Kdr, αr, ηlr, ϕr, kr)
        i = 1

        @views Y    = solid1d_mush_relax.Y[i,:,:]
        @views dYdθ = solid1d_mush_relax.dYdθ[i,:,:]
        @views dYdϕ = solid1d_mush_relax.dYdϕ[i,:,:]
        @views Z    = solid1d_mush_relax.Z[i,:,:]
        @views X    = solid1d_mush_relax.X[i,:,:]

        y1 = y[1]
        y2 = y[2]
        y3 = y[5]
        y4 = y[6]
        y7 = y[4]

        λr = Kdr - 2μr/3
        βr = λr + 2μr

        # Compute strain tensor
        ϵ[:,:,1] = (-2λr*y1 + n*(n+1)λr*y2 + rr*y3 + rr*αr*y7)/(βr*rr)  * Y     # e_rr
        ϵ[:,:,2] = 1/rr * ((y1 - 0.5n*(n+1)y2)Y + 0.5y2*X)                      # e_
        ϵ[:,:,3] = 1/rr * ((y1 - 0.5n*(n+1)y2)Y - 0.5y2*X)                      # e_
        ϵ[:,:,4] = 0.5/μr * y4 * dYdθ                                           # e_rθ
        ϵ[:,:,5] = 0.5/μr * y4 * dYdϕ .* 1.0 ./ sin.(clats)                     # e_rϕ
        ϵ[:,:,6] = 0.5 * y2/rr * Z                                              # e_        
    end


    """
        compute_darcy_displacement!(dis_rel, y, n, r, ω, ϕ, ηl, k, g, ρl)

    Calculate the Darcy displacement vector at a particular radial level. 

    # Arguments
    - `dis_rel::Array{ComplexF64,4}`    : 4D array to store the Darcy displacement vector at a particular radial level, with dimensions corresponding to latitude, longitude, and the 3 components of the Darcy displacement vector.
    - `y::Array{precc,1}`               : 1D array of the solution vector y at a particular radial level, with 8 components.
    - `n::Int`                          : Tidal degree. 
    - `r::prec`                         : Radius at which to compute the Darcy displacement vector.
    - `ω::prec`                         : Forcing frequency.
    - `ϕ::prec`                         : Porosity at radius r.
    - `ηl::prec`                        : Liquid viscosity at radius r.
    - `k::prec`                         : Permeability at radius r.
    - `g::prec`                         : Gravity at radius r.
    - `ρl::prec`                        : Liquid density at radius r.
    """
    function compute_darcy_displacement!(dis_rel, y, n, r, ω, ϕ, ηl, k, g, ρl)
        i = 1

        @views Y    = solid1d_mush_relax.Y[i,:,:]
        @views dYdθ = solid1d_mush_relax.dYdθ[i,:,:]
        @views dYdϕ = solid1d_mush_relax.dYdϕ[i,:,:]
        
        y1 = y[1]
        y5 = y[3]
        y7 = y[4]
        y8 = y[8]
        y9 = 1im * k / (ω*ϕ*ηl*r) * (ρl*g*y1 - ρl * y5 + ρl*g*y8 + y7)
        
        dis_rel[:,:,1] = Y * y8 
        dis_rel[:,:,2] = dYdθ * y9
        dis_rel[:,:,3] = dYdϕ * y9 ./ sin.(clats)
    end


    """
        compute_pore_pressure!(p, y, n)

    Calculate the fluid pore pressur at a particular radial level. 

    # Arguments
    - `p::Array{ComplexF64,4}`          : 4D array to store the pore pressure at a particular radial level, with dimensions corresponding to latitude and longitude.
    - `y::Array{precc,1}`               : 1D array of the solution vector y at a particular radial level, with 8 components.
    - `n::Int`                          : Tidal degree.
    """
    function compute_pore_pressure!(p, y, n)
        i = 1

        @views Y    = solid1d_mush_relax.Y[i,:,:]
        @views dYdθ = solid1d_mush_relax.dYdθ[i,:,:]
        @views dYdϕ = solid1d_mush_relax.dYdϕ[i,:,:]
        
        y7 = y[4]
        
        p[:,:] = Y * y7 
    end


    """
        get_heating_profile(y, r, ρ, g, μ, Ks, ω, ρl, Kl, Kd, α, ηl, ϕ, k, n; lay=nothing)

    Get the radial volumetric heating for two-phase tides and eccentricity forcing,
    assuming synchronous rotation. Heating rate is computed with numerical integration 
    using the solution `y`, using Eq. 2.39a/b/c integrated over solid angle. 

    # Arguments
    - `y::Array{ComplexF64,2}`           : 2D array [state_vector, radius] of the solution vector.
    - `r::AbstractVector`                : 1D vector of radial boundaries/shell coordinates.
    - `ρ::AbstractVector`                : 1D vector of layer densities.
    - `g::AbstractVector`                : 1D vector of gravity values.
    - `μ::AbstractVector`                : 1D vector of complex shear moduli.
    - `Ks::AbstractVector`               : 1D vector of bulk moduli for shear dissipation.
    - `ω::Float64`                       : Tidal frequency in radians per second.
    - `ρl::AbstractVector`               : 1D vector of liquid densities.
    - `Kl::AbstractVector`               : 1D vector of liquid bulk moduli.
    - `Kd::AbstractVector`               : 1D vector of drained bulk moduli.
    - `α::AbstractVector`                : 1D vector of Biot coefficients.
    - `ηl::AbstractVector`               : 1D vector of liquid viscosities.
    - `ϕ::AbstractVector`                : 1D vector of porosities.
    - `k::AbstractVector`                : 1D vector of permeabilities.
    - `n::Int`                           : Tidal degree.

    # Returns
    - `Eμ_tot::Vector{Float64}`          : Total power dissipated due to shear (W/m³).
    - `Eκ_tot::Vector{Float64}`          : Total power dissipated due to compaction (W/m³).
    - `El_tot::Vector{Float64}`          : Total power dissipated due to Darcy flow (W/m³).
    """
    function get_heating_profile(y, r, ρ, g, μ, Ks, ω, ρl, Kl, Kd, α, ηl, ϕ, k, n)

        dres = deg2rad(solid1d_mush_relax.res)
        clats = solid1d_mush_relax.clats
        lons = solid1d_mush_relax.lons
        
        # Reference spherical harmonic Y if needed globally
        # (Assuming the logic uses the same spatial resolution as the first example)
        Y = solid1d_mush_relax.Y[1, :, :] 

        Nr = length(r)
        nlats = length(clats)
        nlons = length(lons)

        # Pre-allocate spatial buffers
        ϵ = zeros(ComplexF64, nlats, nlons, 6)
        d_disp = zeros(ComplexF64, nlats, nlons, 3)
        p = zeros(ComplexF64, nlats, nlons)

        # Output vectors (per radial shell)
        Eμ_tot = zeros(Float64, Nr - 1)
        Eκ_tot = zeros(Float64, Nr - 1)
        El_tot = zeros(Float64, Nr - 1)

        for i in 1:Nr-1
            rr = r[i]
            dr = r[i+1] - r[i]
            dvol = 4π/3 * (r[i+1]^3 - r[i]^3)
            yrr = y[:, i]

            # 1. Compute Tensors/Fields
            compute_strain_ten!(ϵ, yrr, n, rr, ρ[i], g[i], μ[i], Ks[i], ω, ρl[i], Kl[i], Kd[i], α[i], ηl[i], ϕ[i], k[i])
            
            if ϕ[i] > 0
                compute_darcy_displacement!(d_disp, yrr, n, rr, ω, ϕ[i], ηl[i], k[i], g[i], ρl[i])
                # In your snippet, p was being overwritten by y7 * Y immediately
                p .= yrr[7] .* Y 
            else
                p .= 0.0
            end

            # 2. Shear Heating (Solid)
            Eμ_loc = sum(abs.(ϵ[:,:,1:3]).^2, dims=3) .+ 
                     2sum(abs.(ϵ[:,:,4:6]).^2, dims=3) .- 
                     1/3 .* abs.(sum(ϵ[:,:,1:3], dims=3)).^2
            Eμ_loc .*= ω * imag(μ[i])

            # 3. Compaction Heating (Bulk)
            Eκ_loc = ω/2 * imag(Kd[i]) .* abs.(sum(ϵ[:,:,1:3], dims=3)).^2
            if ϕ[i] > 0
                Eκ_loc .+= (ω/2 * imag(Kd[i]) .* (abs.(p) ./ Ks[i]).^2)
            end

            # 4. Darcy Dissipation (Liquid)
            El_loc = zeros(nlats, nlons)
            if ϕ[i] > 0
                El_loc .= 0.5 * ϕ[i]^2 * ω^2 * (ηl[i] / k[i]) .* (abs.(d_disp[:,:,1]).^2 .+ abs.(d_disp[:,:,2]).^2 .+ abs.(d_disp[:,:,3]).^2)
            end

            # 5. Angular Integration and Volume Averaging
            weight = sin.(clats)
            
            Eμ_tot[i] = sum(weight .* Eμ_loc .* dres^2) * rr^2 * dr / dvol
            Eκ_tot[i] = sum(weight .* Eκ_loc .* dres^2) * rr^2 * dr / dvol
            El_tot[i] = sum(weight .* El_loc .* dres^2) * rr^2 * dr / dvol
        end

        return Eμ_tot, Eκ_tot, El_tot
    end

end
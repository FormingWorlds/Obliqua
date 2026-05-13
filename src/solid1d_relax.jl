

module solid1d_relax
    
    include("common.jl")
    using .common

    using LinearAlgebra
    using DoubleFloats
    using AssociatedLegendrePolynomials    
    using StaticArrays
    using SpecialFunctions
    using SparseArrays

    prec  = BigFloat
    precc = Complex{BigFloat}

    # prec  = Float64
    # precc = Complex{Float64}

    const G::prec       = prec(6.6743e-11)       # m^3 kg^-1 s^-2


    """
        resample_profiles(radius, rho, visc, shear, bulk, m_core, dr_min, dr_max)

    Resample the input profiles onto a new grid with `ncalc` points. The new grid is generated using a 
    stretched and refined scheme, which allows for better resolution in regions of interest (e.g., near 
    layer boundaries). 

    # Arguments
    - `radius::Vector{prec}`              : Original radius profile (layer boundaries).
    - `rho::Vector{prec}`                 : Original density profile (defined at layer centers).
    - `visc::Vector{prec}`                : Original viscosity profile (defined at layer centers).
    - `shear::Vector{precc}`              : Original shear modulus profile (defined at layer centers).
    - `bulk::Vector{prec}`                : Original bulk modulus profile (defined at layer centers).
    - `m_core::prec`                      : Mass of the core, used for gravity calculations.
    - `Δr_min::Int64`                     : Minimum grid spacing for the new grid.
    - `Δr_max::Int64`                     : Maximum grid spacing for the new grid.

    # Returns
    Tuple of resampled profiles on the new grid:
    - `r_new_b::Vector{prec}`             : New radius profile at layer boundaries.
    - `ρ_new::Vector{prec}`               : New density profile at layer centers.
    - `η_new::Vector{prec}`               : New viscosity profile at layer centers.
    - `μ_new::Vector{precc}`              : New shear modulus profile at layer centers.
    - `κ_new::Vector{prec}`               : New bulk modulus profile at layer centers.
    - `g_new::Vector{prec}`               : New gravity profile at layer centers.
    - `M_tot::Float64`                    : Total mass enclosed within the outermost layer boundary.
    """ 
    function resample_profiles(radius::Vector{prec}, rho::Vector{prec}, visc::Vector{prec}, shear::Vector{precc}, bulk::Vector{prec}, m_core::prec, dr_min::Int64, dr_max::Int64)
        # setup grids
        α = log(dr_max / dr_min)

        N = Int(round((radius[end] - radius[1]) / dr_min * α / (exp(α) - 1)))

        # indices i = 1:N
        i = collect(1:N)

        # convert to BigFloat for consistency
        i_bf = prec.(i)
        N_bf = prec(N)

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
        κ_new = similar(bulk, N-1)

        for i in 1:N-1
            # find index such that r_b[idx] <= r_new_c[i] < r_b[idx+1]
            idx = searchsortedfirst(radius, r_new_c[i]) - 1
            idx = clamp(idx, 1, length(rho)) # Safety clamp

            ρ_new[i] = rho[idx]
            η_new[i] = visc[idx]
            μ_new[i] = shear[idx]
            κ_new[i] = bulk[idx]
        end

        g_new, M_tot = get_g(r_new_b, ρ_new, m_core) 

        return r_new_b, ρ_new, η_new, μ_new, κ_new, g_new, M_tot
    end


    """
        get_g(r, ρ, m_core)

    Compute the radial gravity structure associated with a density profile `r` at intervals given by `r`.

    # Arguments
    - `r::Array{prec,1}`               : 1D array of layer boundaries. 
    - `ρ::Array{prec,1}`               : 1D array of layer densities. The length of `ρ` must be equal to the number of columns in `r`.
    - `m_core::prec`                   : Mass of the core.

    # Returns
    - `g::Array{prec,1}`               : 1D array of gravity values at the layer boundaries. The dimensions of `g` are the same as `r`.
    - `M_enc::prec`                    : Total mass enclosed within the outermost layer boundary.
    """
    function get_g(r::Vector{prec}, ρ::Vector{prec}, m_core::prec)

        dm = 4.0/3.0 * π .* diff(r.^3) .* ρ

        M_enc = cumsum(dm) .+ m_core
            
        g = G .* M_enc ./ r[2:end].^2

        return g, M_enc[end]
    end


    """
        doublefactorial(n)

    Compute the double factorial of an integer n, defined as n!! = n * (n-2) * (n-4) * ... until 1 or 0.

    # Arguments
    - `n::Integer`                     : The integer for which to compute the double factorial. Must be non-negative.

    # Returns
    - `result::Integer`                 : The double factorial of n.
    """
    function doublefactorial(n::Integer)
        n < 0 && error("doublefactorial not defined for negative n")
        n == 0 && return one(n)
        n == 1 && return one(n)

        result = one(n)
        for k in n:-2:1
            result *= k
        end
        return result
    end


    """
        solve_radial_system(r, ρ, g, μ, K, ω, n, R_planet, ρ_core, μ_core, κ_core, M_tot; core="liquid")

    Solve the radial system of ODEs for the solid-body problem using a relaxation method. This function 
    implements the forward-backward relaxation scheme described in the main text of N. Kobayashi (2006).
    
    # Arguments
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{precc}`                : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `n::Int`                          : Tidal degree.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `μ_core::prec`                    : Shear modulus of the core, used for core boundary conditions.
    - `κ_core::prec`                    : Bulk modulus of the core, used for core boundary conditions.
    - `M_tot::prec`                     : Total mass of the planet, used for gravity calculations.

    # Keyword Arguments
    - `core::String="liquid"            : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.

    # Returns
    - `y_t::Vector{precc}`              : Vector of length 6 representing the tidal solution at the surface (radius = R_planet). This includes the displacements, stresses, and potential at the surface.
    - `y_l::Vector{precc}`              : Vector of length 6 representing the load solution at the surface (radius = R_planet). This includes the displacements, stresses, and potential at the surface.
    - `R::Vector{Matrix{precc}}`        : Vector of 6x6 matrices representing the coefficients of the ODE system at each radial layer.
    - `S::Vector{Matrix{precc}}`        : Vector of 6x6 matrices representing the normalization.
    """    
    function solve_radial_system(r::Vector{prec}, ρ::Vector{prec}, g::Vector{prec}, μ::Vector{precc}, K::Vector{prec}, ω::prec, n::Int, ρ_core::prec, μ_core::prec, κ_core::prec, M_tot::prec; core::String="liquid")

        Nr = length(r)

        # indices for the three components of the relaxation scheme: 
        #   1: core boundary
        #   2: forward propagation
        #   3: surface boundary
        ids = [(1,2), (2, Nr-1), (Nr-1, Nr)]   

        # non-dimensional scaling
        # this implementation needs to be double-checked for consistency.
        R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(prec(1.), prec(1.), prec(1.))
        # R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(r[end], M_tot, g[end])

        ωs = ω * s0
        rs = r ./ R0
        ρs = ρ ./ ρ0
        gs = g ./ g0
        μs = μ ./ μ0
        Ks = K ./ μ0

        R = Vector{Matrix{precc}}(undef, Nr)

        # component 1: apply core boundary condition and get first solution
        C1l, D2l = core_boundary(R, ids[1], rs, ρs, gs, μs, Ks, ωs, ρ_core/ρ0, μ_core/μ0, κ_core/μ0, core, n; G0=G0)

        # component 2: propagate the solution up to the surface (6x6)
        C1l, D2l = propagate_solid(R, C1l, D2l, ids[2], rs, ρs, gs, μs, Ks, ωs, n; G0=G0)

        # component 3: apply surface boundary condition and solve for the final solution at the surface
        y_t, y_l = surface_boundary(R, C1l, D2l, ids[3], rs, ρs, gs, μs, Ks, ωs, n; G0=G0)

        return y_t, y_l, R, S

    end


    """
        core_boundary(R, ids, r, ρ, g, μ, K, ω, ρ_core, μ_core, κ_core, core, n)

    Perform the forward-backward relaxation step at the core boundary. This function implements the recursion described 
    in N. Kobayashi (2007) for the initial step of the relaxation scheme, where we apply the core boundary condition and 
    get the first solution for the first layer above the core.

    # Arguments
    - `R::Vector`                       : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.
    - `r::Vector{prec}`                 : Vector of radial grid points (layer centers).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{precc}                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `μ_core::prec`                    : Shear modulus of the core, used for core boundary conditions.
    - `κ_core::prec`                    : Bulk modulus of the core, used for core boundary conditions.
    - `core::String`                    : Type of core boundary condition to apply.
    - `n::Int`                          : Tidal degree.

    Keyword Arguments
    - `G0::prec=1`                      : Gravitational constant used for non-dimensional scaling.

    # Returns
    - `C1l::Matrix{precc}`              : 3x6 matrix representing the "stored" lower half of the C1 matrix for the next iteration.
    - `D2l::Matrix{precc}`              : 3x6 matrix representing the "stored" lower half of the D2 matrix for the next iteration.
    """
    function core_boundary(R::Vector, ids::Tuple{Int, Int}, r::Vector{prec}, ρ::Vector{prec}, g::Vector{prec}, μ::Vector{precc}, K::Vector{prec}, ω::prec, ρ_core::prec, μ_core::prec, κ_core::prec, core::String, n::Int; G0=1)

        start_id, end_id = ids

        # boundary conditions
        B1 = get_core_bc!(ω, r[start_id], ρ_core, g[start_id], μ_core, κ_core, core, n; G0=G0)        
        
        # first layer (n = 1)
        dr = r[end_id] - r[start_id]

        A1 = get_A(ω, r[start_id], ρ[start_id], g[start_id], μ[start_id], K[start_id], n; G0=G0)
        A2 = get_A(ω, r[end_id], ρ[end_id], g[end_id], μ[end_id], K[end_id], n; G0=G0)

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
        R[start_id] = -S1 \ Q1

        return C1l, D2l
    end

    
    """
        propagate_solid(R, B, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, n)

    Perform the forward-backward relaxation step for the solid propagation segments. This function implements the 
    recursion described in N. Kobayashi (2007) for the segments of the radial grid that correspond to the solid 
    layers, where we propagate the solution up to the surface using the 6x6 system of equations.

    # Arguments
    - `R::Vector`                       : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `B::Vector{Matrix{precc}}`        : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
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

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 3x6 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 3x6 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.
    """
    function propagate_solid(R::Vector, Cn_l::Matrix{precc}, Dnp_l::Matrix{precc}, ids::Tuple{Int, Int}, r::Vector{prec}, ρ::Vector{prec}, g::Vector{prec}, μ::Vector{precc}, K::Vector{prec}, ω::prec, n::Int; G0=1)

        start_id, end_id = ids

        I6 = Matrix{precc}(I, 6, 6)

        Cn_u = zeros(3,6)
        Dnp_u = zeros(3,6)

        # forward recursion
        for i in start_id:end_id

            dr = r[i+1] - r[i]

            # Calculate A at current and next step
            A_n  = get_A(ω, r[i],   ρ[i],   g[i],   μ[i],   K[i],   n; G0=G0)
            A_np = get_A(ω, r[i+1], ρ[i+1], g[i+1], μ[i+1], K[i+1], n; G0=G0)

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
            R[i] = -Xn \ Qn

            # 5. Update the "stored" lower halves for the next iteration
            Cn_l  = Cn[4:6, :]
            Dnp_l = Dnp[4:6, :]
        end

        return Cn_l, Dnp_l
    end


    """
        surface_boundary(R, B, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n, R_planet)

    Perform the forward-backward relaxation step at the surface boundary. This function implements the recursion described 
    in N. Kobayashi (2007) for the final step of the relaxation scheme, where we apply the surface boundary condition and 
    solve for the final solution at the surface, using the 6x6 system of equations.

    # Arguments
    - `R::Vector`                       : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
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
    function surface_boundary(R::Vector, CNm_l::Matrix{precc}, DN_l::Matrix{precc}, ids::Tuple{Int, Int}, r::Vector{prec}, ρ::Vector{prec}, g::Vector{prec}, μ::Vector{precc}, K::Vector{prec}, ω::prec, n::Int; G0=1)

        start_id, end_id = ids

        # tidal surface boundary condition
        BN_t, b_t = get_surface_bc!(r[end], g[end], n, 1, 0, 0, 0; G0=G0)
        # load surface boundary condition
        BN_l, b_l = get_surface_bc!(r[end], g[end], n, 0, 1, 0, 0; G0=G0)

        PN = [CNm_l; zeros(3,6)]
        SN_t = [DN_l; BN_t]
        SN_l = [DN_l; BN_l]

        XN_t = PN * R[start_id] + SN_t
        XN_l = PN * R[start_id] + SN_l

        # solve outer (tides)
        y_t = XN_t \ b_t
        # solve outer (load)
        y_l = XN_l \ b_l

        return y_t, y_l

    end


    """
        get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1, M=6, N=3)

    Get the core boundary condition matrix `B` for the solid-body problem. The core boundary 
    conditions are derived from the requirement that the radial stress at the core-mantle 
    boundary must balance the tidal potential, and that the tangential stresses must vanish.

    # Arguments
    - `ω::prec`                       : Forcing frequency.
    - `r::prec`                          : Radial position of the core-mantle boundary.
    - `ρ::prec`                          : Average core density.
    - `g::prec`                          : Gravity at the core-mantle boundary.
    - `μ::prec`                          : Average core shear modulus.
    - `K::prec`                          : Average core bulk modulus.
    - `type::String`                     : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.
    - `M::Int=6`                         : Number of rows in the Ic matrix. This should be 6 for the solid-body problem.
    - `N::Int=3`                         : Number of linearly independent solutions to compute. This should be 3 for the solid-body problem.

    # Returns
    - `B::Array{precc,2}`                : 3x6 matrix representing the linear relationship between the state variables at the core and the boundary conditions.
    """
    function get_core_bc!(ω::prec, r::prec, ρ::prec, g::prec, μ::prec, K::prec, type::String, n::Int; G0=1)

        Ic = get_Ic(ω, r, ρ, g, μ, K, type, n; G0=G0)

        # Define indices based on Takeuchi & Saito (1972)
        # u vectors (Displacements/Potential): U=1, V=2, phi=5
        # s vectors (Stresses/Potential Flux): X=3, Y=4, psi=6
        idx_u = [1, 2, 5]
        idx_s = [3, 4, 6]

        Mu = Ic[idx_u, :]
        Ms = Ic[idx_s, :]

        # Eq 91
        b = -Mu * pinv(Ms)

        # Construct the 3x6 B matrix
        T = eltype(b)
        B = zeros(T, 3, 6)

        for i in 1:3
            B[i, idx_u[i]] = 1.0               # Identity for u components
            B[i, idx_s[1]] = b[i, 1]           # Coefficient for X
            B[i, idx_s[2]] = b[i, 2]           # Coefficient for Y
            B[i, idx_s[3]] = b[i, 3]           # Coefficient for psi
        end

        return B
    end


    """
        get_surface_bc!(R, g, n, U, U_prime, tau, P; G0=1)

    Get the surface boundary condition vector `b` and matrix `BN` for the solid-body problem. The surface 
    boundary conditions are determined by setting, respectively (U, U', tau, P) to (1,0,0,0) for tidal Love 
    number and (0,1,0,0) for load Love number in system.

    https://hal.science/hal-03421553/document

    # Arguments
    - `R::prec`                          : Mantle radius, used for surface boundary conditions.
    - `g::prec`                          : Gravity at the surface, used for surface boundary conditions.
    - `n::Int`                           : Tidal degree.
    - `U::Int`                           : Tidal potential at the surface.
    - `U_prime::Int`                     : Radial derivative of the tidal potential at the surface.
    - `tau::Int`                         : Tangential tidal stress at the surface.
    - `P::Int`                           : Surface mass load at the surface.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.

    # Returns
    - `BN::Array{precc,2}`               : 3x6 matrix representing the linear relationship between the state variables at the surface and the boundary conditions.
    - `b::Vector{precc}`                 : Vector of length 6 representing the inhomogeneous part of the surface boundary conditions.
    """
    function get_surface_bc!(R::prec, g::prec, n::Int, U::Int, U_prime::Int, tau::Int, P::Int; G0=1)
        
        # Define surface mass load (zeta) based on Farrell/Longman relation
        zeta = ((2 * n + 1) / (4 * pi * G/G0 * R)) * U_prime

        # b vector (Right Hand Side of the B*y = b system)
        b = zeros(precc, 6) 
        
        # radial Stress y3
        b[4] = -g * zeta * (G/G0) / R - P
        
        # tangential Stress y4
        b[5] = tau
        
        # gravitational potential boundary
        b[6] = ((2 * n + 1) / R) * U + 4 * pi * G/G0 * zeta
        
        # construct the 3x6 B matrix
        # this matrix extracts y3, y4, and the combination for y6
        B = zeros(precc, 3, 6)

        # stress components
        B[1, 3] = 1.0  # radial stress y3
        B[2, 4] = 1.0  # tangential stress y4
        
        # potential component
        B[3, 5] = (n + 1) / R
        B[3, 6] = 1.0

        return B, b
    end


    """
        compute_y(r, ρ, g, μ, K, ω, n, R, ρ_core, μ_core, κ_core, M_tot; core="liquid")

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
    - `n::Int`                          : Tidal degree.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.
    - `μ_core::prec`                    : Shear modulus of the core, used for core boundary conditions.
    - `κ_core::prec`                    : Bulk modulus of the core, used for core boundary conditions.
    - `M_tot::prec`                     : Total mass of the planet, used for non-dimensionalization.

    # Keyword Arguments
    - `core::String="liquid"            : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.

    # Returns
    - `y_t::Matrix{ComplexF64}`         : 6xN matrix of the solution at all radial grid points, where N is the number of radial layers. Each column corresponds to a radial grid point, and each row corresponds to a state variable (displacements, stresses, potential).
    - `y_l::Matrix{ComplexF64}`         : 6x1 matrix of the solution at the surface for the load problem.
    """    
    function compute_y(r::Vector{prec}, ρ::Vector{prec}, g::Vector{prec}, μ::Vector{precc}, K::Vector{prec}, ω::prec, n::Int, ρ_core::prec, μ_core::prec, κ_core::prec, M_tot::prec; core="liquid")

        # solve radial system to get surface solution and recursion matrices
        yN_t, yN_l, R, S = solve_radial_system(r, ρ, g, μ, K, ω, n, ρ_core, μ_core, κ_core, M_tot; core=core)

        Nr = length(r)
        T = eltype(yN_t)

        # allocate as 6 x N matrix 
        y_t = Matrix{T}(undef, 6, Nr)
        y_l = Matrix{T}(undef, 6, 1)

        # solve outer  (tides)
        y_t[:, Nr] = yN_t
        # solve outer (load)
        y_l[:, 1] = yN_l

        # back-substitution
        for i in Nr-1:-1:1
            y_t[:, i] = R[i] * y_t[:, i+1]
        end

        # convert to ComplexF64
        y_t = ComplexF64.(y_t)
        y_l = ComplexF64.(y_l)

        return y_t, y_l
    end

end
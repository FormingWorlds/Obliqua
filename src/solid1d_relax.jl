

module solid1d_relax
    
    using LinearAlgebra
    using DoubleFloats
    using AssociatedLegendrePolynomials    
    using StaticArrays
    using SpecialFunctions
    using SparseArrays

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
        resample_profiles(radius, rho, visc, shear, bulk, grav, ncalc)

    Resample the input profiles onto a new grid with `ncalc` points. The new grid is generated using a 
    stretched and refined scheme, which allows for better resolution in regions of interest (e.g., near 
    layer boundaries). 

    # Arguments
    - `radius::Vector{Float64}`           : Original radius profile (layer boundaries).
    - `rho::Vector{Float64}`              : Original density profile (defined at layer centers).
    - `visc::Vector{Float64}`             : Original viscosity profile (defined at layer centers).
    - `shear::Vector{Float64}`            : Original shear modulus profile (defined at layer centers).
    - `bulk::Vector{Float64}`             : Original bulk modulus profile (defined at layer centers).
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
    - `g_new::Vector{prec}`               : New gravity profile at layer centers.
    """ 
    function resample_profiles(radius, rho, visc, shear, bulk, m_core, dr_min, dr_max)
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

        g_new = get_g(r_new_b, ρ_new, m_core) 

        return r_new_b, ρ_new, η_new, μ_new, κ_new, g_new
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
    """
    function get_g(r, ρ, m_core)

        dm = 4.0/3.0 * π .* diff(r.^3) .* ρ

        M_enc = cumsum(dm) .+ m_core
            
        g = G .* M_enc ./ r[2:end].^2

        return g
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
    The grid is internal to solid1d, but can be accessed with 

        solid1d_relax.clats[:] # colatitude grid
        solid1d_relax.lons[:]  # longitude grid
    """
    function define_spherical_grid(res, n, m)
        solid1d_relax.res = res

        # θ and φ grids
        lons = deg2rad.(collect(0:res:360-0.001))'
        clats = deg2rad.(collect(0:res:180))
        clats[1] += 1e-6
        clats[end] -= 1e-6

        # allocate arrays
        solid1d_relax.Y    = zeros(ComplexF64, 1, length(clats), length(lons))
        solid1d_relax.dYdθ = similar(solid1d_relax.Y)
        solid1d_relax.dYdϕ = similar(solid1d_relax.Y)
        solid1d_relax.Z    = similar(solid1d_relax.Y)
        solid1d_relax.X    = similar(solid1d_relax.Y)

        sinθ = sin.(clats)
        cosθ = cos.(clats)
        cotθ = cosθ ./ sinθ
        cscθ = csc.(clats)

        # Normalization factor for spherical harmonics
        norm = sqrt((2*n+1)  * factorial(n-m) / (4π * factorial(n+m)))
        
        i = 1

        # Y
        solid1d_relax.Y[i,:,:] = Ynm(n,m,clats,lons)

        # ∂Y/∂θ
        Pn = Plm.(n, m, cosθ)
        if n > m
            Pn_1 = Plm.(n-1, m, cosθ)
            dPdθ = (n .* cosθ .* Pn .- (n + m) .* Pn_1) ./ (sinθ)
        else
            # m == n -> P_{n-1}^m = 0
            dPdθ = (n .* cosθ .* Pn) ./ (sinθ)
        end
        solid1d_relax.dYdθ[i,:,:] .= dPdθ .* exp.(1im .* m .* lons)

        # ∂Y/∂ϕ
        solid1d_relax.dYdϕ[i,:,:] .= 1im * m .* solid1d_relax.Y[i,:,:]

        # Z = 2 ((1/sinθ) ∂²Y/∂θ∂ϕ - cotθ cscθ ∂Y/∂ϕ)
        solid1d_relax.Z[i,:,:] .= 2 .* (1im * m ./ sinθ .* solid1d_relax.dYdθ[i,:,:] .- cotθ .* cscθ .* solid1d_relax.dYdϕ[i,:,:])

        # X = -2 (cotθ ∂Y/∂θ + csc²θ ∂²Y/∂ϕ²) - n(n+1)) Y
        solid1d_relax.X[i,:,:] .= -2 .* (cotθ .* solid1d_relax.dYdθ[i,:,:] .- cscθ.^2 .* m^2 .* solid1d_relax.Y[i,:,:]) .- n*(n+1) .* solid1d_relax.Y[i,:,:]

        # Normalize
        solid1d_relax.Y[i,:,:]    .*= norm
        solid1d_relax.dYdθ[i,:,:] .*= norm
        solid1d_relax.dYdϕ[i,:,:] .*= norm
        solid1d_relax.Z[i,:,:]    .*= norm
        solid1d_relax.X[i,:,:]    .*= norm

        # save grids
        solid1d_relax.clats = clats
        solid1d_relax.lons  = lons
    end


    """
        get_scales(R0, M0, s0)

    Compute the characteristic scales for the problem based on a reference radius `R0`, mass `M0`, and time 
    scale `s0`. These scales are used to non-dimensionalize the equations and ensure numerical stability.

    # Arguments
    - `R0::prec`                         : Reference radius scale (e.g., planetary radius).
    - `M0::prec`                         : Reference mass scale (e.g., planetary mass).
    - `s0::prec`                         : Reference time scale (e.g., 1 day in seconds).

    # Returns
    Tuple of characteristic scales:
    - `R0::prec`                         : Length scale (m).
    - `M0::prec`                         : Mass scale (kg).
    - `s0::prec`                         : Time scale (s).
    - `ρ0::prec`                         : Density scale (kg/m^3).
    - `G0::prec`                         : Gravitational constant scale (m^3 kg^-1 s^-2).
    - `g0::prec`                         : Gravity scale (m/s^2).
    - `μ0::prec`                         : Shear modulus scale (Pa).
    - `S::Diagonal{prec}`                : Diagonal scaling matrix for state variables.
    - `Sinv::Diagonal{prec}`             : Inverse of the scaling matrix S.
    """
    function get_scales(R0, M0, s0)

        ρ0 = M0 / R0^3            # 1000 kg/m^3
        μ0 = M0 / (R0 * s0^2)     # 1e9 Pa (1 GPa)
        g0 = R0 / s0^2            # 0.1 m/s^2 (Note: check if you want 10 or 0.1)
        G0 = R0^3 / (M0 * s0^2)   # Gravity constant scaling

        S = Diagonal(precc[
            R0,       # y1: radial displacement (m)
            R0,       # y2: tangential displacement (m)
            μ0,       # y3: radial stress (Pa)
            μ0,       # y4: tangential stress (Pa)
            g0*R0,    # y5: potential (m^2/s^2)
            g0        # y6: potential gradient/gravity (m/s^2)
        ])

        Sinv = inv(S)
        return R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv
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
        get_Ic(ω, r, ρ, g, μ, K, type, n; G0=1, M=6, N=3)
            
    Get the core solution vector.
    https://academic.oup.com/gji/article/203/3/2150/2594863
    
    # Arguments
    - `ω::prec`                          : Angular frequency.
    - `r::prec`                          : Radius of the core boundary.
    - `ρ::prec`                          : Density of the core.
    - `g::prec`                          : Gravity at the core boundary.
    - `μ::prec`                          : Shear modulus of the core.
    - `K::prec`                          : Bulk modulus of the core.
    - `type::String`                     : Type of core, either "liquid" or "solid".
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `M::Int=6`                         : Number of rows in the Ic matrix. This should be 6 for the solid-body problem.
    - `N::Int=3`                         : Number of linearly independent solutions to compute. This should be 3 for the solid-body problem.

    # Returns
    - `Ic::Array{precc,2}`               : MxN array of linearly independent solutions at the core boundary. These are used as starting vectors for the numerical integration across the interior.
    """
    function get_Ic(ω, r, ρ, g, μ, K, type, n; G0=1, M=6, N=3)
        Ic = zeros(precc, M, N)

        G_norm = G / G0

        if type=="liquid"
            Ic[1,1] = -r^n / g
            Ic[1,3] = 1.0
            Ic[2,2] = 1.0
            Ic[3,3] = g*ρ
            Ic[5,1] = r^n
            Ic[6,1] = 2(n-1)*r^(n-1)
            Ic[6,3] = 4π * G_norm * ρ 
        elseif type == "inertial"

            @warn "Inertial core boundary conditions have not been fully implemented. Use with caution."

            φ = 4π * G_norm * ρ / 3

            Ic[1,1] = n * r^(n-1)
            Ic[2,1] = r^(n-1)
            Ic[3,1] = 0.0
            Ic[4,1] = 0.0
            Ic[5,1] = -(n*φ - ω^2) * r^n
            Ic[6,1] = -(2*(n-1)*n*φ - (2*n + 1)*ω^2) * r^(n-1)

            α = sqrt(K / ρ)
            f = -ω^2 / φ
            h = f - (n + 1)
            k2 = (ω^2 + 4φ - n*(n+1)*φ^2 / ω^2) / α^2
            k = sqrt(Complex{BigFloat}(k2))
            x = k * r

            x64 = ComplexF64(x)

            jl_n   = sphericalbesselj(n, x64)
            jl_np1 = sphericalbesselj(n+1, x64)

            ϕl   = doublefactorial(2n+1) / x^n * jl_n
            ϕlp1 = doublefactorial(2n+3) / x^(n+1) * jl_np1
            ψl = 2*(2n+3)/x^2 * (1 - ϕl)
            pref = -r^(n+1) / (2n + 3)

            Ic[1,2] = pref * (0.5 * n * h * ψl + f * ϕlp1)
            Ic[2,2] = pref * (0.5 * h * ψl - ϕlp1)
            Ic[3,2] = -φ * r^n * f * ϕl
            Ic[4,2] = 0.0
            Ic[5,2] = -r^(n+2) * (
                (α^2 * f)/r^2 - (3φ*f)/(2*(2n+3)) * ψl
            )
            Ic[6,2] = -r^(n+1) * (
                (2n+1)*(α^2*f)/r^2 -
                (3φ*((2n+1)*f - n*h))/(2*(2n+3)) * ψl
            )

            Ic[:,3] .= 0.0
            Ic[2,3] = 1.0   # tangential slip
        elseif type == "solid" # incompressible solid core
            # First column
            Ic[1, 1] = n*r^( n+1 ) / ( 2*( 2n + 3) )
            Ic[2, 1] = ( n+3 )*r^( n+1 ) / ( 2*( 2n+3 ) * ( n+1 ) )
            Ic[3, 1] = ( n*ρ*g*r + 2*( n^2 - n - 3)*μ ) * r^n / ( 2*( 2n + 3) )
            Ic[4, 1] = n *( n+2 ) * μ * r^n / ( ( 2n + 3 )*( n+1 ) )
            Ic[6, 1] = 2π*G_norm*ρ*n*r^( n+1 ) / ( 2n + 3 )

            # Second column
            Ic[1, 2] = r^( n-1 )
            Ic[2, 2] = r^( n-1 ) / n
            Ic[3, 2] = ( ρ*g*r + 2*( n-1 )*μ ) * r^( n-2 )
            Ic[4, 2] = 2*( n-1 ) * μ * r^( n-2 ) / n
            Ic[6, 2] = 4π*G_norm*ρ*r^( n-1 )

            # Third column
            Ic[3, 3] = -ρ * r^n
            Ic[5, 3] = -r^n
            Ic[6, 3] = -( 2n + 1) * r^( n-1 )
        else
            error("Invalid core type: $type. Must be 'liquid', 'inertial', or 'solid'.")
        end

        return Ic
    end


    """
        get_A(ω, r, ρ, g, μ, K, n; G0=1, λ=nothing)

    Compute the 6x6 `A` matrix in the ODE for the solid-body problem.

    # Arguments
    - `ω::prec`                          : Forcing frequency of the tidal forcing.
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.
    - `λ::prec=nothing`                  : Lamé's first parameter at radius r. If not provided, it is computed as λ = K - 2μ/3.

    # Returns
    - `A::Array{precc,2}`               : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.
    """
    function get_A(ω, r, ρ, g, μ, K, n; G0=1, λ=nothing)
        A = zeros(precc, 6, 6) 
        get_A!(A, ω, r, ρ, g, μ, K, n; G0=G0, λ=λ)
        return A
    end


    """
        get_A!(A, ω, r, ρ, g, μ, K, n; G0=1, λ=nothing)

    Compute the 6x6 `A` matrix in the ODE for the solid-body problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025) when α=φ=0, as well as Sabadini and Vermeersen 
    (2016) Eq. 1.95.

    # Arguments
    - `A::Array{precc,2}`                : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.
    - `ω::prec`                          : Forcing frequency of the tidal forcing.
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.
    - `λ::prec=nothing`                  : Lamé's first parameter at radius r. If not provided, it is computed as λ = K - 2μ/3.
    """
    function get_A!(A::Matrix, ω, r, ρ, g, μ, K, n; G0=1, λ=nothing)
        if isnothing(λ)
            λ = K - 2μ/3
        end

        G_norm = G / G0

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)
        rβ_inv = r_inv * β_inv

        A[1,1] = -2λ * r_inv*β_inv
        A[2,1] = -r_inv
        A[3,1] = 4r_inv * (3K*μ*r_inv*β_inv - ρ*g) - ω^2 * ρ 
        A[4,1] = -r_inv * (6K*μ*r_inv*β_inv - ρ*g )
        A[5,1] = 4π * G_norm * ρ
        A[6,1] = 4π*(n+1)*G_norm*ρ*r_inv

        A[1,2] = n*(n+1) * λ * r_inv*β_inv
        A[2,2] = r_inv
        A[3,2] = -n*(n+1)*r_inv * (6K*μ*r_inv*β_inv - ρ*g ) 
        A[4,2] = 2μ*r_inv^2 * (n*(n+1)*(1 + λ*β_inv) - 1.0 ) - ω^2 * ρ 
        A[6,2] = -4π*n*(n+1)*G_norm*ρ*r_inv

        A[1,3] = β_inv
        A[3,3] = r_inv*β_inv * (-4μ )
        A[4,3] = -λ * r_inv*β_inv
        
        A[2,4] = 1.0 / μ
        A[3,4] = n*(n+1)*r_inv
        A[4,4] = -3r_inv

        A[3,5] = ρ * (n+1)*r_inv
        A[4,5] = -ρ*r_inv
        A[5,5] = -(n+1)r_inv     

        A[3,6] = -ρ
        A[5,6] = 1.0
        A[6,6] = (n-1)r_inv
    end


    """
        solve_radial_system(r, ρ, g, μ, K, ω, n, R_planet, ρ_core; core="liquid")

    Solve the radial system of ODEs for the solid-body problem using a relaxation method. This function 
    implements the forward-backward relaxation scheme described in the main text of N. Kobayashi (2006).
    
    # Arguments
    - `r::Vector{prec}`                 : Vector of radial grid points (layer boundaries).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `n::Int`                          : Tidal degree.
    - `R_planet::prec`                  : Planetary radius, used for surface boundary conditions.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.

    # Keyword Arguments
    - `core::String="liquid"            : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.

    # Returns
    - `XN::Matrix{precc}`               : 6x6 matrix representing the solution at the surface (radius = R_planet). Each column corresponds to a linearly independent solution of the ODE system, and each row corresponds to a state variable (displacements, stresses, potential).
    - `R::Vector{Matrix{precc}}`        : Vector of 6x6 matrices representing the solution at each radial grid point. Each matrix contains the state variables (displacements, stresses, potential) for the 6-component system of ODEs.
    - `b::Vector{precc}`                : Vector of length 6 representing the inhomogeneous part of the surface boundary conditions.
    """
    function solve_radial_system(r, ρ, g, μ, K, ω, n, R_planet, ρ_core; core="liquid")

        Nr = length(r)

        # non-dimensional scaling
        # this implementation needs to be double-checked for consistency.
        R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(9.8e6, 4.3e22, 4.3e3)

        ωs = ω * s0
        rs = r ./ R0
        ρs = ρ ./ ρ0
        gs = g ./ g0
        μs = μ ./ μ0
        Ks = K ./ μ0

        # boundary conditions
        B1    = get_core_bc!(ω, r[1], ρ_core, g[1], μ[1], K[1], core, n; G0=1)        
        BN, b = get_surface_bc!(R_planet, n)

        # storage
        R = Vector{Matrix{precc}}(undef, Nr)

        # first layer (n = 1)
        dr = rs[2] - rs[1]

        A1 = S * get_A(ωs, rs[1], ρs[1], gs[1], μs[1], Ks[1], n; G0=G0) * Sinv
        A2 = S * get_A(ωs, rs[2], ρs[2], gs[2], μs[2], Ks[2], n; G0=G0) * Sinv

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
        R[1] = -pinv(S1) * Q1

        # forward recursion
        for i in 2:Nr-1

            dr = rs[i+1] - rs[i]

            A_n  = S * get_A(ωs, rs[i],   ρs[i],   gs[i],   μs[i],   Ks[i],   n; G0=G0) * Sinv
            A_np = S * get_A(ωs, rs[i+1], ρs[i+1], gs[i+1], μs[i+1], Ks[i+1], n; G0=G0) * Sinv

            Cn  =  I6 + 0.5 * dr * A_n
            Dnp = -I6 + 0.5 * dr * A_np

            # split
            Cn_u, Cn_l = Cn[1:3, :], Cn[4:6, :]
            Dnp_u, Dnp_l = Dnp[1:3, :], Dnp[4:6, :]

            # build blocks
            Pn = [Cn_l; zeros(3,6)]
            Sn = [Dnp_l; Cn_u]
            Qn = [zeros(3,6); Dnp_u]

            # recursion
            Xn = Pn * R[i-1] + Sn
            R[i] = -pinv(Xn) * Qn
        end

        # final layer (n = N)
        dr = rs[Nr] - rs[Nr-1]

        A_Nm = S * get_A(ωs, rs[end-1], ρs[end-1], gs[end-1], μs[end-1], Ks[end-1], n; G0=G0) * Sinv
        A_N  = S * get_A(ωs, rs[end],   ρs[end],   gs[end],   μs[end],   Ks[end],   n; G0=G0) * Sinv

        CNm =  I6 + 0.5 * dr * A_Nm
        DN  = -I6 + 0.5 * dr * A_N

        CNm_l = CNm[4:6, :]
        DN_l  = DN[4:6, :]

        PN = [CNm_l; zeros(3,6)]
        SN = [DN_l; BN]

        XN = PN * R[Nr-1] + SN

        return XN, R, b
    end


    """
        get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1)

    Get the core boundary condition matrix `B` and vector `b` for the solid-body problem. The core 
    boundary conditions are derived from the requirement that the radial stress at the core-mantle 
    boundary must balance the tidal potential, and that the tangential stresses must vanish.

    # Arguments
    - `ω::Float64`                       : Forcing frequency.
    - `r::prec`                          : Radial position of the core-mantle boundary.
    - `ρ::prec`                          : Density at the core-mantle boundary.
    - `g::prec`                          : Gravity at the core-mantle boundary.
    - `μ::precc`                         : Complex shear modulus at the core-mantle boundary.
    - `K::prec`                          : Bulk modulus at the core-mantle boundary.
    - `type::String`                     : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.

    # Returns
    - `B::Array{precc,2}`                : 3x6 matrix representing the linear relationship between the state variables at the core and the boundary conditions.
    - `b::Vector{precc}`                 : Vector of length 6 representing the inhomogeneous part of the core boundary conditions.
    """
    function get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1)

        Ic = get_Ic(ω, r, ρ, g, μ, K, type, n; G0=G0, M=6, N=3)

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
        get_surface_bc!(R, n)

    Get the surface boundary condition vector `b` and matrix `BN` for the solid-body problem. The surface 
    boundary conditions are derived from the requirement that the radial stress at the surface must balance 
    the tidal potential, and that the tangential stresses must vanish.

    # Arguments
    - `R::prec`                          : Planetary radius, used for surface boundary conditions.
    - `n::Int`                           : Tidal degree.

    # Returns
    - `BN::Array{precc,2}`               : 3x6 matrix representing the linear relationship between the state variables at the surface and the boundary conditions.
    - `b::Vector{precc}`                 : Vector of length 6 representing the inhomogeneous part of the surface boundary conditions.
    """
    function get_surface_bc!(R, n)
        
        # first three elements of b are zero, last three 
        # elements are given in the documentation.
        b = zeros(precc, 6) 
        b[6] = (2n+1)/R
        
        # s vectors (Stresses/Potential Flux): X=3, Y=4, psi=6
        idx_s = [3, 4, 6]

        # Construct the 3x6 B matrix
        B = zeros(3, 6)

        for i in 1:3
            B[i, idx_s[i]] = 1.0               # Identity for s components
        end

        return B, b
    end


    """
        compute_y_relaxation(r, ρ, g, μ, K, ω, n, R, ρ_core; core="liquid")

    Compute the solution `y` to the solid-body problem using a relaxation method. This function performs the 
    forward-backward relaxation scheme described in the main text of N. Kobayashi (2006), where we first solve 
    the radial system of ODEs to obtain the solution at the surface, and then perform back-substitution to 
    compute the solution at all interior points. 

    # Arguments
    - `r::Vector{prec}`                 : Vector of radial grid points (layer boundaries).
    - `ρ::Vector{prec}`                 : Vector of densities at the layer centers.
    - `g::Vector{prec}`                 : Vector of gravity values at the layer centers.
    - `μ::Vector{prec}`                 : Vector of shear moduli at the layer centers.
    - `K::Vector{prec}`                 : Vector of bulk moduli at the layer centers.
    - `ω::prec`                         : Angular frequency of the tidal forcing.
    - `n::Int`                          : Tidal degree.
    - `R::prec`                         : Planetary radius, used for surface boundary conditions.
    - `ρ_core::prec`                    : Density of the core, used for core boundary conditions.

    # Keyword Arguments
    - `core::String="liquid"            : Type of core boundary condition to apply. Options are "liquid" for a fluid core, "solid" for a solid core, and "inertial" for a core with inertial response.

    # Returns
    - `y::Matrix{precc}`                : 6xN matrix of the solution at all radial grid points, where N is the number of radial layers. Each column corresponds to a radial grid point, and each row corresponds to a state variable (displacements, stresses, potential).
    """    
    function compute_y_relaxation(r, ρ, g, μ, K, ω, n, R, ρ_core; core="liquid")

        # solve radial system to get surface solution and recursion matrices
        XN, R, b = solve_radial_system(r, ρ, g, μ, K, ω, n, R, ρ_core; core=core)

        Nr = length(r)
        T = eltype(XN)

        # allocate as 6 x N matrix 
        y = Matrix{T}(undef, 6, Nr)

        # solve outer boundary
        y[:, Nr] = pinv(XN) * b

        # back-substitution
        for i in Nr-1:-1:1
            y[:, i] = R[i] * y[:, i+1]
        end

        return y
    end


    """
        compute_strain_ten!(ϵ, y, n, rr, ρr, gr, μr, Ksr)

    Calculate the strain tensor ϵ at a particular radial level. 

    # Arguments
    - `ϵ::Array{ComplexF64,3}`          : 3D array to store the strain tensor at a particular radial level, with dimensions corresponding to latitude, longitude, and the 6 independent components of the strain tensor.
    - `y::Array{precc,1}`               : 1D array of the solution vector y at a particular radial level, with 6 components.
    - `n::Int`                          : Tidal degree.
    - `rr::prec`                        : Radius at which to compute the strain tensor.
    - `ρr::prec`                        : Density at radius rr.
    - `gr::prec`                        : Gravity at radius rr.
    - `μr::prec`                        : Shear modulus at radius rr.
    - `Ksr::prec`                       : Bulk modulus at radius rr.
    """
    function compute_strain_ten!(ϵ, y, n, rr, ρr, gr, μr, Ksr)
               
        i = 1

        @views Y    = solid1d_relax.Y[i,:,:]
        @views dYdθ = solid1d_relax.dYdθ[i,:,:]
        @views dYdϕ = solid1d_relax.dYdϕ[i,:,:]
        @views Z    = solid1d_relax.Z[i,:,:]
        @views X    = solid1d_relax.X[i,:,:]

        y1 = y[1]
        y2 = y[2]
        y3 = y[3]
        y4 = y[4]

        λr = Ksr .- 2μr/3
        βr = λr + 2μr

        # Compute strain tensor
        ϵ[:,:,1] = (-2λr*y1 + n*(n+1)λr*y2 + rr*y3)/(βr*rr)  * Y
        ϵ[:,:,2] = 1/rr * ((y1 - 0.5n*(n+1)y2)Y + 0.5y2*X)
        ϵ[:,:,3] = 1/rr * ((y1 - 0.5n*(n+1)y2)Y - 0.5y2*X)
        ϵ[:,:,4] = 0.5/μr * y4 * dYdθ        
        ϵ[:,:,5] = 0.5/μr * y4 * dYdϕ .* 1.0 ./ sin.(clats) 
        ϵ[:,:,6] = 0.5 * y2/rr * Z
    end


    """
        function get_heating_profile(y, r, ρ, g, μ, κ, n, ω; lay=nothing)

    Get the radial volumetric heating for solid-body tides and eccentricity forcing,
    assuming synchronous rotation. Heating rate is computed with numerical integration 
    using the solution `y` returned by [`compute_y`](@ref), using Eq. 2.39a/b integrated 
    over solid angle. The heating profile for a specific layer is specified with `lay`, 
    otherwise all layers will be caclulated.

    # Arguments
    - `y::Array{ComplexF64,4}`           : 4D array of the solution vector y across the interior, returned by `compute_y`.
    - `r::Array{Float64,2}`              : 2D array of layer boundaries.
    - `ρ::Array{Float64,1}`              : 1D array of layer densities.
    - `g::Array{Float64,2}`              : 2D array of gravity values at the layer boundaries.
    - `μ::Array{Float64,1}`              : 1D array of layer shear moduli.
    - `κ::Array{Float64,1}`              : 1D array of layer bulk moduli.
    - `n::Int`                           : Tidal degree.
    - `ω::Float64`                       : Tidal frequency in radians per second.

    # Keyword Arguments
    - `lay::Int=nothing`                 : If specified, compute the heating profile for only the layer corresponding to this index. Otherwise, compute for all layers.

    # Returns
    - `Eμ_tot::Array{Float64,1}`         : 1D array of total power dissipated in each primary layer due to shear, in W.
    - `Eκ_tot::Array{Float64,1}`         : 1D array of total power dissipated in each primary layer due to compaction, in W.
    """
    function get_heating_profile(y, r, ρ, g, μ, κ, n, ω)

        dres = deg2rad(solid1d_relax.res)

        clats = solid1d_relax.clats
        lons  = solid1d_relax.lons

        Nr = length(r)
        nlats = length(clats)
        nlons = length(lons)

        # strain tensor per radius
        ϵ = zeros(ComplexF64, nlats, nlons, 6)

        # outputs
        Eμ_vol = zeros(Float64, Nr-1)
        Eκ_vol = zeros(Float64, Nr-1)

        Eμ_tot = zeros(Float64, Nr-1)
        Eκ_tot = zeros(Float64, Nr-1)

        for i in 1:Nr-1

            rr = r[i]
            dr = r[i+1] - r[i]

            dvol = 4π/3 * (r[i+1]^3 - r[i]^3)

            yrr = y[:, i]

            compute_strain_ten!(ϵ, yrr, n, rr, ρ[i], g[i], μ[i], κ[i])

            # shear heating
            Eμ_loc = sum(abs.(ϵ[:,:,1:3]).^2, dims=3) .+
                    2sum(abs.(ϵ[:,:,4:6]).^2, dims=3) .-
                    1/3 .* abs.(sum(ϵ[:,:,1:3], dims=3)).^2

            Eμ_loc .*= ω * imag(μ[i])

            # bulk heating
            Eκ_loc = ω/2 * imag(κ[i]) .* abs.(sum(ϵ[:,:,1:3], dims=3)).^2

            # angular integration
            weight = sin.(clats)

            Eμ_tot[i] = sum(weight .* Eμ_loc * dres^2) * rr^2 * dr / dvol
            Eκ_tot[i] = sum(weight .* Eκ_loc * dres^2) * rr^2 * dr / dvol

        end

        return Eμ_tot, Eκ_tot
    end

end
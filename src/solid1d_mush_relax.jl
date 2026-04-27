

module solid1d_mush_relax
    
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
        get_scales(R0, M0, g0)

    Compute the characteristic scales for the problem based on a reference radius `R0`, mass `M0`, and gravity scale 
    `g0`. These scales are used to non-dimensionalize the equations and ensure numerical stability.

    # Arguments
    - `R0::prec`                         : Reference radius scale (e.g., planetary radius).
    - `M0::prec`                         : Reference mass scale (e.g., planetary mass).
    - `g0::prec`                         : Reference gravity scale.

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
    function get_scales(R0, M0, g0)

        ρ0 = M0 / R0^3
        P0 = g0 * R0
        μ0 = ρ0 * g0 * R0

        s0 = sqrt(g0 / R0)
        G0 = R0^3 / (M0 * s0^2) 

        S = Diagonal(precc[
            R0,       # y1: radial displacement (m)
            R0,       # y2: tangential displacement (m)
            P0,    # y5: potential (m^2/s^2)
            μ0,       # y7: pore pressure (Pa)
            μ0,       # y3: radial stress (Pa)
            μ0,       # y4: tangential stress (Pa)
            g0,       # y6: potential gradient/gravity (m/s^2)
            R0,       # y8: relative radial displacement (m)
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
        get_Ic(ω, r, ρ, g, μ, K, type, n; G0=1, M=8, N=4)
            
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
    function get_Ic(ω, r, ρ, g, μ, K, type, n; G0=1, M=8, N=4)
        Ic = zeros(precc, M, N)

        G_norm = G / G0

        if type=="liquid"
            if M == 6
                Ic[1,1] = -r^n / g
                Ic[1,3] = 1.0
                Ic[2,2] = 1.0
                Ic[4,3] = g*ρ
                Ic[3,1] = r^n
                Ic[6,1] = 2(n-1)*r^(n-1)
                Ic[6,3] = 4π * G_norm * ρ 
            elseif M == 8
                Ic[1,1] = -r^n / g
                Ic[1,3] = 1.0
                Ic[2,2] = 1.0
                Ic[5,3] = g*ρ
                Ic[3,1] = r^n
                Ic[7,1] = 2(n-1)*r^(n-1)
                Ic[7,3] = 4π * G_norm * ρ 
            end
        elseif type == "inertial"
            @warn "Inertial core boundary conditions have not been fully implemented. Use with caution."
            
            φ = 4π * G_norm * ρ / 3
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

            if M == 6
                Ic[1,1] = n * r^(n-1)
                Ic[2,1] = r^(n-1)
                Ic[3,1] = 0.0
                Ic[4,1] = 0.0
                Ic[5,1] = -(n*φ - ω^2) * r^n
                Ic[6,1] = -(2*(n-1)*n*φ - (2*n + 1)*ω^2) * r^(n-1)

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
            elseif M == 8
                Ic[1,1] = n * r^(n-1)
                Ic[2,1] = r^(n-1)
                Ic[5,1] = 0.0
                Ic[6,1] = 0.0
                Ic[3,1] = -(n*φ - ω^2) * r^n
                Ic[7,1] = -(2*(n-1)*n*φ - (2*n + 1)*ω^2) * r^(n-1)

                Ic[1,2] = pref * (0.5 * n * h * ψl + f * ϕlp1)
                Ic[2,2] = pref * (0.5 * h * ψl - ϕlp1)
                Ic[5,2] = -φ * r^n * f * ϕl
                Ic[6,2] = 0.0
                Ic[3,2] = -r^(n+2) * (
                    (α^2 * f)/r^2 - (3φ*f)/(2*(2n+3)) * ψl
                )
                Ic[7,2] = -r^(n+1) * (
                    (2n+1)*(α^2*f)/r^2 -
                    (3φ*((2n+1)*f - n*h))/(2*(2n+3)) * ψl
                )

                Ic[:,3] .= 0.0
                Ic[2,3] = 1.0   # tangential slip
            end
        elseif type == "solid" # incompressible solid core
            if M == 6
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
            elseif M == 8
                # First column
                Ic[1, 1] = n*r^( n+1 ) / ( 2*( 2n + 3) )
                Ic[2, 1] = ( n+3 )*r^( n+1 ) / ( 2*( 2n+3 ) * ( n+1 ) )
                Ic[5, 1] = ( n*ρ*g*r + 2*( n^2 - n - 3)*μ ) * r^n / ( 2*( 2n + 3) )
                Ic[6, 1] = n *( n+2 ) * μ * r^n / ( ( 2n + 3 )*( n+1 ) )
                Ic[7, 1] = 2π*G_norm*ρ*n*r^( n+1 ) / ( 2n + 3 )

                # Second column
                Ic[1, 2] = r^( n-1 )
                Ic[2, 2] = r^( n-1 ) / n
                Ic[5, 2] = ( ρ*g*r + 2*( n-1 )*μ ) * r^( n-2 )
                Ic[6, 2] = 2*( n-1 ) * μ * r^( n-2 ) / n
                Ic[7, 2] = 4π*G_norm*ρ*r^( n-1 )

                # Third column
                Ic[5, 3] = -ρ * r^n
                Ic[3, 3] = -r^n
                Ic[7, 3] = -( 2n + 1) * r^( n-1 )
            end
        else
            error("Invalid core type: $type. Must be 'liquid', 'inertial', or 'solid'.")
        end

        return Ic
    end


    """
        get_A(ω, r, ρ, g, μ, K, n; G0=1, λ=nothing, M=8)

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
    - `M::Int=8`                         : Number of rows in the A matrix. This should be 6 for the solid-body problem, but can be 8 for the two-phase problem.

    # Returns
    - `A::Array{precc,2}`               : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.
    """
    function get_A(ω, r, ρ, g, μ, K, n; G0=1, λ=nothing, M=8)
        A = zeros(precc, 6, 6) 
        get_A!(A, ω, r, ρ, g, μ, K, n; G0=G0, λ=λ, M=M)
        return A
    end

    
    """
        get_A(ω, r, ρ, g, μ, K, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1, λ=nothing)

    Compute the 8x8 `A` matrix in the ODE for the two-phase problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025).

    # Arguments
    - `ω::prec`                          : Forcing frequency.
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `ρₗ::prec`                          : Liquid density at radius r.
    - `Kl::prec`                         : Liquid bulk modulus at radius r.
    - `Kd::prec`                         : Drained bulk modulus at radius r.
    - `α::prec`                          : Biot coefficient at radius r.
    - `ηₗ::prec`                          : Liquid viscosity at radius r.
    - `ϕ::prec`                          : Porosity at radius r.
    - `k::prec`                          : Permeability at radius r.    
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `G0::prec=1`                       : Gravitational constant scale for non-dimensionalization.
    - `λ::prec=nothing`                  : Lamé's first parameter at radius r. If not provided, it is computed as λ = K - 2μ/3.

    # Returns
    - `A::Array{precc,2}`               : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.

    See also [`get_A!`](@ref)
    """
    function get_A(ω, r, ρ, g, μ, K, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1, λ=nothing)
        A = zeros(precc, 8, 8)
        get_A!(A, ω, r, ρ, g, μ, K, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=G0, λ=λ)
        return A
    end


    """
        get_A!(A, ω, r, ρ, g, μ, K, n; G0=1, λ=nothing, M=8)

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
    - `M::Int=8`                         : Number of rows in the A matrix. This should be 6 for the solid-body problem, but can be 8 for the two-phase problem.
    """
    function get_A!(A::Matrix, ω, r, ρ, g, μ, K, n; G0=1, λ=nothing, M=8)
        if isnothing(λ)
            λ = K - 2μ/3
        end

        G_norm = G / G0

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)
        rβ_inv = r_inv * β_inv

        if M == 8
            A[1,1] = -2λ * r_inv*β_inv
            A[2,1] = -r_inv
            A[5,1] = 4r_inv * (3K*μ*r_inv*β_inv - ρ*g) - ω^2 * ρ 
            A[6,1] = -r_inv * (6K*μ*r_inv*β_inv - ρ*g )
            A[3,1] = 4π * G_norm * ρ
            A[7,1] = 4π*(n+1)*G_norm*ρ*r_inv

            A[1,2] = n*(n+1) * λ * r_inv*β_inv
            A[2,2] = r_inv
            A[5,2] = -n*(n+1)*r_inv * (6K*μ*r_inv*β_inv - ρ*g ) 
            A[6,2] = 2μ*r_inv^2 * (n*(n+1)*(1 + λ*β_inv) - 1.0 ) - ω^2 * ρ 
            A[7,2] = -4π*n*(n+1)*G_norm*ρ*r_inv

            A[1,5] = β_inv
            A[5,5] = r_inv*β_inv * (-4μ )
            A[6,5] = -λ * r_inv*β_inv
            
            A[2,6] = 1.0 / μ
            A[5,6] = n*(n+1)*r_inv
            A[6,6] = -3r_inv

            A[5,3] = ρ * (n+1)*r_inv
            A[6,3] = -ρ*r_inv
            A[3,3] = -(n+1)r_inv     

            A[5,7] = -ρ
            A[3,7] = 1.0
            A[7,7] = (n-1)r_inv

        elseif M ==6
            A[1,1] = -2λ * r_inv*β_inv
            A[2,1] = -r_inv
            A[4,1] = 4r_inv * (3K*μ*r_inv*β_inv - ρ*g) - ω^2 * ρ 
            A[5,1] = -r_inv * (6K*μ*r_inv*β_inv - ρ*g )
            A[3,1] = 4π * G_norm * ρ
            A[6,1] = 4π*(n+1)*G_norm*ρ*r_inv

            A[1,2] = n*(n+1) * λ * r_inv*β_inv
            A[2,2] = r_inv
            A[4,2] = -n*(n+1)*r_inv * (6K*μ*r_inv*β_inv - ρ*g ) 
            A[5,2] = 2μ*r_inv^2 * (n*(n+1)*(1 + λ*β_inv) - 1.0 ) - ω^2 * ρ 
            A[6,2] = -4π*n*(n+1)*G_norm*ρ*r_inv

            A[1,4] = β_inv
            A[4,4] = r_inv*β_inv * (-4μ )
            A[5,4] = -λ * r_inv*β_inv
            
            A[2,5] = 1.0 / μ
            A[4,5] = n*(n+1)*r_inv
            A[5,5] = -3r_inv

            A[4,3] = ρ * (n+1)*r_inv
            A[5,3] = -ρ*r_inv
            A[3,3] = -(n+1)r_inv     

            A[4,6] = -ρ
            A[3,6] = 1.0
            A[6,6] = (n-1)r_inv
        end
    end


    """
        get_A!(A, ω, r, ρ, g, μ, K, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Compute the 8x8 `A` matrix in the ODE for the two-phase problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025).

    # Arguments
    - `A::Array{precc,2}`                : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.
    - `ω::prec`                          : Forcing frequency.
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `ρₗ::prec`                          : Liquid density at radius r.
    - `Kl::prec`                         : Liquid bulk modulus at radius r.
    - `Kd::prec`                         : Drained bulk modulus at radius r.
    - `α::prec`                          : Biot coefficient at radius r.
    - `ηₗ::prec`                          : Liquid viscosity at radius r.
    - `ϕ::prec`                          : Porosity at radius r.
    - `k::prec`                          : Permeability at radius r.
    - `n::Int`                           : Tidal degree.

    # Notes
    See also [`get_A`](@ref)
    """
    function get_A!(A::Matrix, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1, λ=nothing)
        # λ = K - 2μ/3       # Lame's second param, which uses the drained compaction modulus
        λ = Kd .- 2μ/3       # Lame's second param, which uses the drained compaction modulus
        S = ϕ/Kl + (α - ϕ)/K # Storavity, which uses liquid and solid grain bulk moduli  

        # First add the solid-body coefficients, but using drained moduli. 
        get_A!(A, ω, r, ρ, g, μ, Kd, n; λ=λ, G0=G0, M=8)    # Note that here we replace the bulk modulus with the compaction modulus

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)

        G_norm = G / G0

        # ϕ = 0.
        # If there is a porous layer, now add the two-phase components
        if !iszero(ϕ)

            A[1,4] = α * β_inv

            A[5,1] += 1im * k*ρₗ^2 *g^2 * n*(n+1) / (ω*ηₗ) * r_inv^2
            A[5,3] += -(n+1)r_inv * 1im *(k*ρₗ^2*g*n)/(ω*ηₗ) * r_inv                               
            A[5,4] = 1im * (k*ρₗ*g*n*(n+1))/(ω*ηₗ)*r_inv^2 - 4μ*α*β_inv*r_inv
            A[5,8] =  1im * (k*ρₗ^2*g^2*n*(n+1))/(ω*ηₗ)*r_inv^2 - 4ϕ*ρₗ*g*r_inv 
        
            A[6,4] = 2α*μ*r_inv * β_inv
            A[6,8] = ϕ*ρₗ*g*r_inv 
            
            A[3,8] = 4π*G_norm*ρₗ*ϕ

            A[7,1] += -1im * 4π*G_norm*n*(n+1)*r_inv * (k*ρₗ^2*g/(ω*ηₗ)*r_inv)
            A[7,3] = 1im*4π*n*(n+1)G_norm*(ρₗ)^2*k*r_inv^2 / (ω*ηₗ)
            A[7,4] = -1im *4π*n*(n+1)G_norm*ρₗ*k*r_inv^2 / ( ω*ηₗ) 
            A[7,8] = 4π*G_norm*(n+1)*r_inv * (ϕ*ρₗ - 1im * n*k*ρₗ^2*g/(ω*ηₗ)*r_inv) 
            
            A[4,1] = ρₗ*g*r_inv * ( 4 - 1im *(k*ρₗ*g*n*(n+1)/(ω*ϕ*ηₗ))*r_inv)  
            A[4,2] = -ρₗ*n*(n+1)*r_inv*g
            A[4,3] = -ρₗ*(n+1)r_inv * (1 - 1im*(k*ρₗ*g*n)/(ω*ϕ*ηₗ)*r_inv)  
            A[4,7] = ρₗ 
            A[4,4] = - 1im*(k*ρₗ*g*n*(n+1))/(ω*ϕ*ηₗ)*r_inv^2
            A[4,8] = -1im*ω*ϕ*ηₗ/k - 4π*G_norm*(ρ - ϕ*ρₗ)*ρₗ + ρₗ*g*r_inv*(4 - 1im*(k*ρₗ*g*n*(n+1))/(ω*ϕ*ηₗ)*r_inv) 
        
            A[8,1] = r_inv*( 1im * k*ρₗ*g*n*(n+1)/(ω*ϕ*ηₗ)*r_inv - α/ϕ * 4μ*β_inv) 
            A[8,2] = α/ϕ * 2n*(n+1)*μ *β_inv * r_inv
            A[8,5] = -α/ϕ * β_inv 
            A[8,3] = -1im * k *ρₗ *n*(n+1) / (ω*ϕ*ηₗ)*r_inv^2 
            A[8,4] = 1im*k*n*(n+1)/(ω*ϕ*ηₗ)*r_inv^2 - 1/ϕ * (S + α^2 * β_inv) # If solid and liquid are compressible, keep the 1/M term
            A[8,8] = 1im * k *ρₗ*g *n*(n+1) / (ω*ϕ*ηₗ)*r_inv^2  - 2r_inv 
        end

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
    - `B::Vector{Matrix{precc}}`        : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
    - `S::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the normalization.
    - `ifc::Int`                        : Index of the first interface layer (the one closer to the core).
    - `ifd::Int`                        : Index of the second interface layer (the one closer to the surface).
    """
    function solve_radial_system(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n,
                                ρ_core, M_tot; core="liquid")

        # 1. Find the original interface index
        ifc_orig = findfirst(k .> 0)
        ifd_orig = findlast(k .> 0)

        vars = (r, ρ, g, μ, K, ρₗ, Kl, Kd, α, ηₗ, ϕ, k)
        new_vars = map(vars) do v
            v_new = copy(v)
            insert!(v_new, ifc_orig, v[ifc_orig])
            insert!(v_new, ifd_orig + 1, v[ifd_orig])
            v_new
        end
        new_r, new_ρ, new_g, new_μ, new_K, new_ρₗ, new_Kl, new_Kd, new_α, new_ηₗ, new_ϕ, new_k = new_vars

        Nr = length(new_r)
        ifc = ifc_orig  # The first of the two duplicate layers
        ifd = ifd_orig + 1  # The second of the two duplicate layers

        # 3. Define the new segments for the relaxation scheme
        # Segment 3 now covers the transition between the two identical radial points
        ids = [
            (1, 2),             # 1: Core Boundary
            (2, ifc-1),         # 2: Solid Propagation
            (ifc-1, ifc+1),     # 3: Interface Transition (duplicate layer)
            (ifc+1, ifd-1),     # 4: Mushy Propagation
            (ifd-1, ifd+1),     # 5: Interface Transition (duplicate layer)
            (ifd+1, Nr-1),      # 6: Solid Propagation
            (Nr-1, Nr)          # 7: Surface Boundary
        ]

        # 4. Non-dimensional scaling (Not working atm, use 1.,1.,1. for now)
        R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(1., 1., 1.)
        # R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(r[end], M_tot, g[end])
        
        ωs = ω * s0
        rs = new_r ./ R0
        ρs = new_ρ ./ ρ0
        gs = new_g ./ g0
        μs = new_μ ./ μ0
        Ks = new_K ./ μ0
        ρₗs = new_ρₗ ./ ρ0
        Kls = new_Kl ./ μ0
        Kds = new_Kd ./ μ0
        ηₗs = new_ηₗ ./ (μ0 * s0)
        ks = new_k ./ R0^2

        # 5. Initialize Matrices
        R = [zeros(precc, 8, 8) for _ in 1:Nr]
        B = [zeros(precc, 8, 1) for _ in 1:Nr]

        # Define the specific indices used by 6x6
        idx = [1, 2, 3, 5, 6, 7]

        # Create the view using these indices for both rows and columns
        R6_view = [view(R[i], idx, idx) for i in 1:Nr]
        R8_view = [view(R[i], 1:8, 1:8) for i in 1:Nr]
        B6_view = [view(B[i], idx, 1) for i in 1:Nr]
        B8_view = [view(B[i], 1:8, 1) for i in 1:Nr]

        # component 1: apply core boundary condition and get first solution (3x6)
        C1l, D2l = core_boundary(R6_view, ids[1], rs, ρs, gs, μs, Ks, ωs, ρ_core/ρ0, core, n; G0=G0)

        # component 1: apply core boundary condition and get first solution (4x8)
        # C1l, D2l = core_boundary_mush(R8_view, ids[1], rs, ρs, gs, μs, Ks, ωs, ρₗs, Kls, Kds, α, ηₗs, ϕ, ks, ρ_core/ρ0, core, n; G0=G0)

        # component 2: propagate the solution up to the surface (6x6)
        C1l, D2l = propagate_solid(R6_view, B6_view, C1l, D2l, ids[2], rs, ρs, gs, μs, Ks, ωs, n; G0=G0)

        # component 3: interface between 6x6 and 8x8
        C1l, D2l = interface_solid_mush(R8_view, B8_view, C1l, D2l, ids[3])

        # component 4: propagate the solution up to the surface (8x8)
        C1l, D2l = propagate_mush(R8_view, B8_view, C1l, D2l, ids[4], rs, ρs, gs, μs, Ks, ωs, ρₗs, Kls, Kds, new_α, ηₗs, new_ϕ, ks, n; G0=G0)

        # component 5: interface between 8x8 and 6x6
        C1l, D2l = interface_mush_solid(R8_view, B8_view, C1l, D2l, ids[5])

        # component 6: propagate the solution up to the surface (6x6)
        C1l, D2l = propagate_solid(R6_view, B6_view, C1l, D2l, ids[6], rs, ρs, gs, μs, Ks, ωs, n; G0=G0)

        # component 7: apply surface boundary condition and solve for the final solution at the surface
        y_t, y_l = surface_boundary(R6_view, B6_view, C1l, D2l, ids[7], rs, ρs, gs, μs, Ks, ωs, n; G0=G0)

        # component 3: apply surface boundary condition and solve for the final solution at the surface
        # y_t, y_l = surface_boundary_mush(R8_view, B8_view, C1l, D2l, ids[5], rs, ρs, gs, μs, Ks, ωs, n; G0=G0)

        return y_t, y_l, R, B, S, ifc, ifd
    end


    """
        interface_mush_solid(R8, B8, Cn_l, Dnp_l, ids)

    Perform the forward-backward relaxation step at the interface between the mushy layer and the solid layer. This 
    function implements the recursion described in N. Kobayashi (2007) for the transition from the 8x8 system to the 6x6 system.

    # Arguments
    - `R8::Vector{Matrix{precc}}`       : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `B8::Vector{Matrix{precc}}`       : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
    - `Cn_l::Matrix{precc}`             : 3x6 matrix representing the "stored" lower half of the Cn matrix from the previous step.
    - `Dnp_l::Matrix{precc}`            : 3x6 matrix representing the "stored" lower half of the Dnp matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 3x6 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 3x6 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.
    """
    function interface_mush_solid(R8, B8, Cn_l, Dnp_l, ids)

        start_id, end_id = ids

        # Impose continuity at the boundary
        # Make sure bn[4, 1] = 0.0, i.e. zero darcy flux at boundary
        bn = zeros(precc, 8, 1) 

        I86 = zeros(precc, 8, 8)
        I8  = Matrix{precc}(I, 8, 8)

        icc = [1.,1.,1.,0.,1.,1.,1.,1.]
        idd = [1.,1.,1.,0.,1.,1.,1.,0.]
        for i in 1:8
            I86[i, i] = icc[i]
            I8[i, i]  = idd[i]
        end      

        Cn  = I86
        Dnp = -I8

        target_cols = [1, 2, 3, 5, 6, 7]
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

        # 4. Perform recursion
        Xn = Pn * R8[start_id] + Sn
        R_ifc = - pinv(Xn) * Qn

        R8[start_id+1] .= R_ifc 
        B8[start_id+1] .= pinv(Xn) * (bn - Pn * B8[start_id])

        # 5. Update the "stored" lower halves for the next iteration
        Cn_l  = Cn[5:7, target_cols]
        Dnp_l = Dnp[5:7, target_cols]

        return Cn_l, Dnp_l
    end


    """
        interface_solid_mush(R8, B8, Cn_l, Dnp_l, ids)

    Perform the forward-backward relaxation step at the interface between the solid layer and the mushy layer. This 
    function implements the recursion described in N. Kobayashi (2007) for the transition from the 6x6 system to the 8x8 system.

    # Arguments
    - `R8::Vector{Matrix{precc}}`       : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `B8::Vector{Matrix{precc}}`       : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
    - `Cn_l::Matrix{precc}`             : 3x6 matrix representing the "stored" lower half of the Cn matrix from the previous step.
    - `Dnp_l::Matrix{precc}`            : 3x6 matrix representing the "stored" lower half of the Dnp matrix from the previous step.
    - `ids::Tuple{Int, Int}`            : Tuple containing the start and end indices of the current segment in the radial grid.

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 3x6 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 3x6 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.    
    """
    function interface_solid_mush(R8, B8, Cn_l, Dnp_l, ids)

        start_id, end_id = ids

        # impose continuity at the boundary
        bn = zeros(precc, 8, 1) 
        # Induce some pore pressure 
        bn[4, 1] = -1.0

        I86 = zeros(precc, 8, 8)
        iss = [1.,1.,1.,1.,1.,1.,1.,0.]
        for i in 1:8
            I86[i, i] = iss[i]
        end

        I8  = Matrix{precc}(I, 8, 8)

        Cn  = I86
        Dnp = -I8

        target_cols = [1, 2, 3, 5, 6, 7]
        # 1. Use the "stored" lower halves from the previous step 
        # to fill the upper blocks of P and S.
        Pn_u = zeros(precc, 4, 8)
        Pn_u[1:3, target_cols] .= Cn_l

        Sn_u = zeros(precc, 4, 8)
        Sn_u[1:3, target_cols] .= Dnp_l

        Qn_u = zeros(precc, 4, 8)

        # 2. Get the upper halves of the NEWLY calculated Cn and Dnp
        Cn_u  = Cn[1:4, :]
        Dnp_u = Dnp[1:4, :]

        # 3. Build the 8x8 blocks
        Pn = [Pn_u; zeros(precc, 4, 8)]
        Sn = [Sn_u; Cn_u]
        Qn = [Qn_u; Dnp_u]

        # 4. Perform recursion
        Xn = Pn * R8[start_id] + Sn
        R_ifc = - pinv(Xn) * Qn

        R8[start_id+1] .= R_ifc 
        B8[start_id+1] .= pinv(Xn) * (bn - Pn * B8[start_id])

        # 5. Update the "stored" lower halves for the next iteration
        Cn_l  = Cn[5:8, :]
        Dnp_l = Dnp[5:8, :]

        return Cn_l, Dnp_l

    end


    """
        core_boundary(R, ids, r, ρ, g, μ, K, ω, ρ_core, core, n)

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

    # Returns
    - `C1l::Matrix{precc}`              : 3x6 matrix representing the "stored" lower half of the C1 matrix for the next iteration.
    - `D2l::Matrix{precc}`              : 3x6 matrix representing the "stored" lower half of the D2 matrix for the next iteration.
    """
    function core_boundary(R, ids, r, ρ, g, μ, K, ω, ρ_core, core, n; G0=1)

        start_id, end_id = ids

        # boundary conditions
        B1 = get_core_bc!(ω, r[start_id], ρ_core, g[start_id], μ[start_id], K[start_id], core, n; G0=G0, M=6, N=3)        
        
        # first layer (n = 1)
        dr = r[end_id] - r[start_id]

        A1 = get_A(ω, r[start_id], ρ[start_id], g[start_id], μ[start_id], K[start_id], n; G0=G0, M=6)
        A2 = get_A(ω, r[end_id], ρ[end_id], g[end_id], μ[end_id], K[end_id], n; G0=G0, M=6)

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
        core_boundary_mush(R, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, ρ_core, core, n)

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
    function core_boundary_mush(R, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, ρ_core, core, n; G0=1)

        start_id, end_id = ids

        # boundary conditions
        B1 = get_core_bc!(ω, r[start_id], ρ_core, g[start_id], μ[start_id], K[start_id], core, n; G0=G0, M=8, N=4)     
        
        # first layer (n = 1)
        dr = r[end_id] - r[start_id]

        A1 = get_A(ω, r[start_id], ρ[start_id], g[start_id], μ[start_id], K[start_id],
                    ρₗ[start_id], Kl[start_id], Kd[start_id], α[start_id], ηₗ[start_id], ϕ[start_id], k[start_id], n; G0=G0)

        A2 = get_A(ω, r[end_id], ρ[end_id], g[end_id], μ[end_id], K[end_id],
                ρₗ[end_id], Kl[end_id], Kd[end_id], α[end_id], ηₗ[end_id], ϕ[end_id], k[end_id], n; G0=G0)

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
        propagate_solid(R, B, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, n)

    Perform the forward-backward relaxation step for the solid propagation segments. This function implements the 
    recursion described in N. Kobayashi (2007) for the segments of the radial grid that correspond to the solid 
    layers, where we propagate the solution up to the surface using the 6x6 system of equations.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
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
    function propagate_solid(R, B, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, n; G0=1)

        start_id, end_id = ids

        I6 = Matrix{precc}(I, 6, 6)

        Cn_u = zeros(3,6)
        Dnp_u = zeros(3,6)

        # forward recursion
        for i in start_id:end_id

            dr = r[i+1] - r[i]

            # Calculate A at current and next step
            A_n  = get_A(ω, r[i],   ρ[i],   g[i],   μ[i],   K[i],   n; G0=G0, M=6)
            A_np = get_A(ω, r[i+1], ρ[i+1], g[i+1], μ[i+1], K[i+1], n; G0=G0, M=6)

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
            R[i] .= -Xn \ Qn
            B[i] .=  Xn \ (-Pn * B[i-1])

            # 5. Update the "stored" lower halves for the next iteration
            Cn_l  = Cn[4:6, :]
            Dnp_l = Dnp[4:6, :]
        end

        return Cn_l, Dnp_l
    end

    
    """
        propagate_mush(R, B, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Perform the forward-backward relaxation step for the mushy layer propagation segment. This function implements the
    recursion described in N. Kobayashi (2007) for the segment of the radial grid that corresponds to the mushy layer, 
    where we propagate the solution up to the surface using the full 8x8 system of equations that includes the porous 
    layer components.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `B::Vector{Matrix{precc}}`        : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
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

    # Returns
    - `Cn_l::Matrix{precc}`             : Updated 4x8 matrix representing the "stored" lower half of the Cn matrix for the next iteration.
    - `Dnp_l::Matrix{precc}`            : Updated 4x8 matrix representing the "stored" lower half of the Dnp matrix for the next iteration.
    """
    function propagate_mush(R, B, Cn_l, Dnp_l, ids, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; G0=1)

        start_id, end_id = ids

        I8 = Matrix{precc}(I, 8, 8)

        Cn_u = zeros(4,8)
        Dnp_u = zeros(4,8)

        # forward recursion
        for i in start_id:end_id

            dr = r[i+1] - r[i]

            # Calculate A at current and next step
            A_n  = get_A(ω, r[i],   ρ[i],   g[i],   μ[i],   K[i],
                            ρₗ[i], Kl[i], Kd[i], α[i], ηₗ[i], ϕ[i], k[i], n; G0=G0)

            A_np = get_A(ω, r[i+1], ρ[i+1], g[i+1], μ[i+1], K[i+1],
                        ρₗ[i+1], Kl[i+1], Kd[i+1], α[i+1], ηₗ[i+1], ϕ[i+1], k[i+1], n; G0=G0)

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

            R[i] .= -Xn \ Qn
            B[i] .=  Xn \ (-Pn * B[i-1])

            # 5. Update the "stored" lower halves for the next iteration
            Cn_l  = Cn[5:8, :]
            Dnp_l = Dnp[5:8, :]
        end

        return Cn_l, Dnp_l
    end


    """
        surface_boundary(R, B, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n)

    Perform the forward-backward relaxation step at the surface boundary. This function implements the recursion described 
    in N. Kobayashi (2007) for the final step of the relaxation scheme, where we apply the surface boundary condition and 
    solve for the final solution at the surface, using the 6x6 system of equations.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `B::Vector{Matrix{precc}}`        : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
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
    function surface_boundary(R, B, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n; G0=1)

        start_id, end_id = ids

        # tidal surface boundary condition
        BN_t, b_t = get_surface_bc!(r[end], g[end], n, 1, 0, 0, 0; M=6, N=3)
        # load surface boundary condition
        BN_l, b_l = get_surface_bc!(r[end], g[end], n, 0, 1, 0, 0; M=6, N=3)

        PN = [CNm_l; zeros(3,6)]
        SN_t = [DN_l; BN_t]
        SN_l = [DN_l; BN_l]

        XN_t = PN * R[start_id] + SN_t
        XN_l = PN * R[start_id] + SN_l

        BN = - XN_t \ PN * B[start_id]

        # solve outer  (tides)
        y_t = XN_t \ b_t + BN
        # solve outer (load)
        y_l = XN_l \ b_l + BN

        return y_t, y_l

    end


    """
        surface_boundary_mush(R, B, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n)

    Perform the forward-backward relaxation step at the surface boundary for the two-phase problem. This function implements 
    the recursion described in N. Kobayashi (2007) for the final step of the relaxation scheme, where we apply the surface 
    boundary condition and solve for the final solution at the surface, using the full 8x8 system of equations that includes 
    the porous layer components.

    # Arguments
    - `R::Vector{Matrix{precc}}`        : Vector of 8x8 matrices representing the coefficients of the ODE system at each radial layer.
    - `B::Vector{Matrix{precc}}`        : Vector of 8x1 matrices representing the inhomogeneous terms of the ODE system at each radial layer.
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
    function surface_boundary_mush(R, B, CNm_l, DN_l, ids, r, ρ, g, μ, K, ω, n; G0=1)

        start_id, end_id = ids

        # tidal surface boundary condition
        BN_t, b_t = get_surface_bc!(r[end], g[end], n, 1, 0, 0, 0; M=8, N=4)
        # load surface boundary condition
        BN_l, b_l = get_surface_bc!(r[end], g[end], n, 0, 1, 0, 0; M=8, N=4)

        PN = [CNm_l; zeros(4,8)]
        SN_t = [DN_l; BN_t]
        SN_l = [DN_l; BN_l]

        XN_t = PN * R[start_id] + SN_t
        XN_l = PN * R[start_id] + SN_l

        BN = - XN_t \ PN * B[start_id]

        # solve outer  (tides)
        y_t = pinv(XN_t) * b_t + BN
        # solve outer (load)
        y_l = pinv(XN_l) * b_l + BN

        return y_t, y_l

    end


    """
        get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1, M=6, N=3)

    Get the core boundary condition matrix `B` for the solid-body problem. The core boundary 
    conditions are derived from the requirement that the radial stress at the core-mantle 
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
    - `M::Int=6`                         : Dimensionality of the system (6 for standard, 8 for mushy layer).
    - `N::Int=3`                         : Number of boundary conditions to apply (3 for standard, 4 for mushy layer).

    # Returns
    - `B::Array{precc,2}`                : 3x6 matrix representing the linear relationship between the state variables at the core and the boundary conditions.
    """
    function get_core_bc!(ω, r, ρ, g, μ, K, type, n; G0=1, M=6, N=3)
        # 1. Get the Initial Conditions matrix
        Ic = get_Ic(ω, r, ρ, g, μ, K, type, n; G0=G0, M=M, N=N)

        # 2. Define indices based on dimensionality
        # If M=8 (Mushing/Hay 2025):  U=1, V=2, phi=5, P=7 | X=3, Y=4, psi=6, R=8
        # If M=6 (Standard/Takeuchi): U=1, V=2, phi=5      | X=3, Y=4, psi=6
        if M == 8
            idx_u = [1, 2, 3, 4]
            idx_s = [5, 6, 7, 8]
        else
            idx_u = [1, 2, 5]
            idx_s = [3, 4, 6]
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
        get_surface_bc!(R, g, n, U, U_prime, tau, P; G0=1, M=6, N=3)

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
    - `M::Int=6`                         : Dimensionality of the system (6 for standard, 8 for mushy layer).
    - `N::Int=3`                         : Number of boundary conditions to apply (3 for standard, 4 for mushy layer).

    # Returns
    - `B::Array{precc,2}`                : 3x6 matrix representing the linear relationship between the state variables at the surface and the boundary conditions.
    - `b::Vector{precc}`                 : Vector of length 6 representing the inhomogeneous part of the surface boundary conditions.
    """
    function get_surface_bc!(R, g, n, U, U_prime, tau, P; G0=1, M=8, N=4)
        
        # Define surface mass load (zeta) based on Farrell/Longman relation
        zeta = ((2 * n + 1) / (4 * pi * G/G0 * R)) * U_prime

        # b vector (Right Hand Side of the B*y = b system)
        b = zeros(precc, M) 
        
        if M == 8
            # radial Stress y3
            b[5] = -g * zeta * (G/G0) / R - P
            
            # tangential Stress y4
            b[6] = tau
            
            # gravitational potential boundary
            b[7] = ((2 * n + 1) / R) * U + 4 * pi * G/G0 * zeta

            # darcy flux boundary
            b[8] = 0.
        elseif M == 6
            # radial Stress y3
            b[4] = -g * zeta * (G/G0) / R - P
            
            # tangential Stress y4
            b[5] = tau
            
            # gravitational potential boundary
            b[6] = ((2 * n + 1) / R) * U - 4 * pi * G/G0 * zeta
        else
            error("Unsupported M value. M should be either 6 or 8.")
        end
        
        # construct the 4x8 B matrix
        # this matrix extracts y3, y4, and the combination for y6
        B = zeros(precc, N, M)

        if M == 8
            # stress components
            B[1, 5] = 1.0  # radial stress y3
            B[2, 6] = 1.0  # tangential stress y4
            
            # potential component
            B[3, 3] = (n + 1) / R
            B[3, 7] = 1.0
            B[4, 8] = 1.0
        elseif M == 6
            # stress components
            B[1, 4] = 1.0  # radial stress y3
            B[2, 5] = 1.0  # tangential stress y4
            # potential component
            B[3, 3] = (n + 1) / R
            B[3, 6] = 1.0        
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
        yN_t, yN_l, R, B, S, ifc, ifd = solve_radial_system(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n, ρ_core, M_tot; core=core)

        Nr = length(r) + 2
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
            if i > ifc && i < ifd
                # in the mushy region, use the mushy recursion
                y_t[:, i] = R[i] * y_t[:, i+1] + B[i]
            else
                # in the solid region, use the solid recursion
                y_t[:, i] = R[i] * y_t[:, i+1] + B[i]
            end
        end

        # Create a mask for all columns except the two interface indices
        mask = [i for i in 1:size(y_t, 2) if i != ifc && i != ifd]

        # Apply the mask to the columns
        y_t = y_t[:, mask]

        # apply scaling to get dimensional solution
        for i in 1:Nr-2
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


module solid1d_mush
    
    using LinearAlgebra
    using DoubleFloats
    using AssociatedLegendrePolynomials    
    using StaticArrays


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
        expand_layers(r; nr::Int=80)

    Discretize the primary layers given by `r` into `nr` discrete secondary layers.

    # Arguments
    - `r::Array{Float64,2}`               : 2D array of primary layer boundaries.

    # Keyword Arguments
    - `nr::Int=80`                        : Number of secondary layers to discretize.

    # Returns
    - `rs::Array{Float64,2}`              : 2D array of secondary layer boundaries/
    """
    function expand_layers(r; nr::Int=80)
        
        rs = zeros(Float64, (nr+1, length(r)-1))
        
        for i in 1:length(r)-1
            rfine = LinRange(r[i], r[i+1], nr+1)
            rs[:, i] .= rfine[1:end] 
        end
    
        return rs
    end


    """
        get_g(r, ρ)

    Compute the radial gravity structure associated with a density profile `r` at intervals given by `r`.

    # Arguments
    - `r::Array{Float64,2}`               : 2D array of layer boundaries. 
    - `ρ::Array{Float64,1}`               : 1D array of layer densities. The length of `ρ` must be equal to the number of columns in `r`.

    # Returns
    - `g::Array{Float64,2}`               : 2D array of gravity values at the layer boundaries. The dimensions of `g` are the same as `r`.

    # Notes
    `r` must be be a 2D array, with index 1 representing the top radius of secondary layers, and index 2
    representing the top radius of primary layers. 
    """
    function get_g(r, ρ)
        g = zeros(Float64, size(r))
        M = zeros(Float64, size(r))

        for i in 1:size(r)[2]
            M[2:end,i] = 4.0/3.0 * π .* diff(r[:,i].^3) .* ρ[i]
        end
    
        g[2:end,:] .= G*accumulate(+,M[2:end,:]) ./ r[2:end,:].^2
        g[1,2:end] = g[end,1:end-1]

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

        solid1d.clats[:] # colatitude grid
        solid1d.lons[:]  # longitude grid
    """
    function define_spherical_grid(res, n, m)
        solid1d_mush.res = res

        # θ and φ grids
        lons = deg2rad.(collect(0:res:360-0.001))'
        clats = deg2rad.(collect(0:res:180))
        clats[1] += 1e-6
        clats[end] -= 1e-6

        # allocate arrays
        solid1d_mush.Y    = zeros(ComplexF64, 1, length(clats), length(lons))
        solid1d_mush.dYdθ = similar(solid1d_mush.Y)
        solid1d_mush.dYdϕ = similar(solid1d_mush.Y)
        solid1d_mush.Z    = similar(solid1d_mush.Y)
        solid1d_mush.X    = similar(solid1d_mush.Y)

        sinθ = sin.(clats)
        cosθ = cos.(clats)
        cotθ = cosθ ./ sinθ
        cscθ = csc.(clats)

        # Normalization factor for spherical harmonics
        norm = sqrt((2*n+1)  * factorial(n-m) / (4π * factorial(n+m)))
        
        i = 1

        # Y
        solid1d_mush.Y[i,:,:] = Ynm(n,m,clats,lons)

        # ∂Y/∂θ
        Pn = Plm.(n, m, cosθ)
        if n > m
            Pn_1 = Plm.(n-1, m, cosθ)
            dPdθ = (n .* cosθ .* Pn .- (n + m) .* Pn_1) ./ (sinθ)
        else
            # m == n -> P_{n-1}^m = 0
            dPdθ = (n .* cosθ .* Pn) ./ (sinθ)
        end
        solid1d_mush.dYdθ[i,:,:] .= dPdθ .* exp.(1im .* m .* lons)

        # ∂Y/∂ϕ
        solid1d_mush.dYdϕ[i,:,:] .= 1im * m .* solid1d_mush.Y[i,:,:]

        # Z = 2 ((1/sinθ) ∂²Y/∂θ∂ϕ - cotθ cscθ ∂Y/∂ϕ)
        solid1d_mush.Z[i,:,:] .= 2 .* (1im * m ./ sinθ .* solid1d_mush.dYdθ[i,:,:] .- cotθ .* cscθ .* solid1d_mush.dYdϕ[i,:,:])

        # X = -2 (cotθ ∂Y/∂θ + csc²θ ∂²Y/∂ϕ²) - n(n+1)) Y
        solid1d_mush.X[i,:,:] .= -2 .* (cotθ .* solid1d_mush.dYdθ[i,:,:] .- cscθ.^2 .* m^2 .* solid1d_mush.Y[i,:,:]) .- n*(n+1) .* solid1d_mush.Y[i,:,:]

        # Normalize
        solid1d_mush.Y[i,:,:]    .*= norm
        solid1d_mush.dYdθ[i,:,:] .*= norm
        solid1d_mush.dYdϕ[i,:,:] .*= norm
        solid1d_mush.Z[i,:,:]    .*= norm
        solid1d_mush.X[i,:,:]    .*= norm

        # save grids
        solid1d_mush.clats = clats
        solid1d_mush.lons  = lons
    end


    """
        convert_params_to_prec(r, ρ, g, μ, κs, ω, ρl, κl, κd, α, ηl, ϕ, k)

    Convert input parameters into the required precision.
    # Arguments
    - `r::Array{Float64,2}`               : 2D array of layer boundaries.
    - `ρ::Array{Float64,1}`               : 1D array of layer densities. 
    - `g::Array{Float64,2}`               : 2D array of gravity values at the layer boundaries. 
    - `μ::Array{Float64,1}`               : 1D array of layer shear moduli.
    - `κs::Array{Float64,1}`              : 1D array of layer bulk moduli. 
    - `ω::Float64`                        : Forcing frequency.
    - `ρl::Array{Float64,1}`              : 1D array of layer liquid densities.
    - `κl::Array{Float64,1}`              : 1D array of layer liquid bulk moduli.
    - `κd::Array{Float64,1}`              : 1D array of layer drained bulk moduli.
    - `α::Array{Float64,1}`               : 1D array of layer Biot coefficients.
    - `ηl::Array{Float64,1}`              : 1D array of layer liquid viscosities.
    - `ϕ::Array{Float64,1}`               : 1D array of layer porosities.
    - `k::Array{Float64,1}`               : 1D array of layer permeabilities.

    # Returns
    Tuple of converted parameters in the required precision.
    """
    function convert_params_to_prec(r, ρ, g, μ, κs, ω, ρl, κl, κd, α, ηl, ϕ, k)
        r_prec = convert(Array{prec}, r)
        ρ_prec = convert(Array{prec}, ρ)
        g_prec = convert(Array{prec}, g)
        μ_prec = convert(Array{precc}, μ)
        κs_prec = convert(Array{precc}, κs)
        ρl_prec = convert(Array{prec}, ρl)
        κl_prec = convert(Array{prec}, κl)
        κd_prec = convert(Array{precc}, κd)
        α_prec = convert(Array{precc}, α)
        ηl_prec = convert(Array{prec}, ηl)
        ϕ_prec = convert(Array{prec}, ϕ)
        k_prec = convert(Array{prec}, k)
        ω_prec = convert(prec, ω)

        return (r_prec,  ρ_prec, g_prec, μ_prec, κs_prec, ω_prec, ρl_prec, κl_prec, κd_prec, α_prec, ηl_prec, ϕ_prec, k_prec)
    end


    """
        get_Ic(r, ρ, g, μ, type, n; M=8, N=4)
            
    Get the core solution vector.
    
    # Arguments
    - `r::prec`                          : Radius of the core boundary.
    - `ρ::prec`                          : Density of the core.
    - `g::prec`                          : Gravity at the core boundary.
    - `μ::prec`                          : Shear modulus of the core.
    - `type::String`                     : Type of core, either "liquid" or "solid".
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `M::Int=8`                         : Number of rows in the Ic matrix. This should be 8 for the two-phase problem.
    - `N::Int=4`                         : Number of linearly independent solutions to compute. This should be 4 for the two-phase problem.

    # Returns
    - `Ic::Array{precc,2}`               : MxN array of linearly independent solutions at the core boundary. These are used as starting vectors for the numerical integration across the interior.
    """
    function get_Ic(r, ρ, g, μ, type, n; M=8, N=4)
        Ic = zeros(precc, M, N)

        if type=="liquid"
            Ic[1,1] = -r^n / g
            Ic[1,3] = 1.0
            Ic[2,2] = 1.0
            Ic[3,3] = g*ρ
            Ic[5,1] = r^n
            Ic[6,1] = 2(n-1)*r^(n-1)
            Ic[6,3] = 4π * G * ρ 
        else # incompressible solid core
            # First column
            Ic[1, 1] = n*r^( n+1 ) / ( 2*( 2n + 3) )
            Ic[2, 1] = ( n+3 )*r^( n+1 ) / ( 2*( 2n+3 ) * ( n+1 ) )
            Ic[3, 1] = ( n*ρ*g*r + 2*( n^2 - n - 3)*μ ) * r^n / ( 2*( 2n + 3) )
            Ic[4, 1] = n *( n+2 ) * μ * r^n / ( ( 2n + 3 )*( n+1 ) )
            Ic[6, 1] = 2π*G*ρ*n*r^( n+1 ) / ( 2n + 3 )

            # Second column
            Ic[1, 2] = r^( n-1 )
            Ic[2, 2] = r^( n-1 ) / n
            Ic[3, 2] = ( ρ*g*r + 2*( n-1 )*μ ) * r^( n-2 )
            Ic[4, 2] = 2*( n-1 ) * μ * r^( n-2 ) / n
            Ic[6, 2] = 4π*G*ρ*r^( n-1 )

            # Third column
            Ic[3, 3] = -ρ * r^n
            Ic[5, 3] = -r^n
            Ic[6, 3] = -( 2n + 1) * r^( n-1 )

        end

        return Ic
    end


    """
        get_A(r, ρ, g, μ, K, n)

    Compute the 6x6 `A` matrix in the ODE for the solid-body problem.

    # Arguments
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `n::Int`                           : Tidal degree.

    # Returns
    - `A::Array{precc,2}`               : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.

    # Notes
    See also [`get_A!`](@ref)
    """
    function get_A(r, ρ, g, μ, K, n)
        A = zeros(precc, 6, 6) 
        get_A!(A, r, ρ, g, μ, K, n)
        return A
    end


    """
        get_A(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Compute the 8x8 `A` matrix in the ODE for the two-phase problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025).

    # Arguments
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `ω::prec`                          : Forcing frequency.
    - `ρₗ::prec`                          : Liquid density at radius r.
    - `Kl::prec`                         : Liquid bulk modulus at radius r.
    - `Kd::prec`                         : Drained bulk modulus at radius r.
    - `α::prec`                          : Biot coefficient at radius r.
    - `ηₗ::prec`                          : Liquid viscosity at radius r.
    - `ϕ::prec`                          : Porosity at radius r.
    - `k::prec`                          : Permeability at radius r.    
    - `n::Int`                           : Tidal degree.

    # Returns
    - `A::Array{precc,2}`               : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.

    See also [`get_A!`](@ref)
    """
    function get_A(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        A = zeros(precc, 8, 8)
        get_A!(A, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        return A
    end


    """
        get_A!(A, r, ρ, g, μ, K, n; λ=nothing)

    Compute the 6x6 `A` matrix in the ODE for the solid-body problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025) when α=φ=0, as well as Sabadini and Vermeersen 
    (2016) Eq. 1.95.

    # Arguments
    - `A::Array{precc,2}`                : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `λ::prec=nothing`                  : Lamé's first parameter at radius r. If not provided, it is computed as λ = K - 2μ/3.

    # Notes
    See also [`get_A`](@ref)
    """
    function get_A!(A::Matrix, r, ρ, g, μ, K, n; λ=nothing)
        if isnothing(λ)
            λ = K - 2μ/3
        end

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)
        rβ_inv = r_inv * β_inv

        A[1,1] = -2λ * r_inv*β_inv
        A[2,1] = -r_inv
        A[3,1] = 4r_inv * (3K*μ*r_inv*β_inv - ρ*g)       #- ω^2 * ρ# 
        A[4,1] = -r_inv * (6K*μ*r_inv*β_inv - ρ*g )
        A[5,1] = 4π * G * ρ
        A[6,1] = 4π*(n+1)*G*ρ*r_inv

        A[1,2] = n*(n+1) * λ * r_inv*β_inv
        A[2,2] = r_inv
        A[3,2] = -n*(n+1)*r_inv * (6K*μ*r_inv*β_inv - ρ*g ) 
        A[4,2] = 2μ*r_inv^2 * (n*(n+1)*(1 + λ*β_inv) - 1.0 ) #- ω^2 * ρ# 
        A[6,2] = -4π*n*(n+1)*G*ρ*r_inv

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
        get_A!(A, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Compute the 8x8 `A` matrix in the ODE for the two-phase problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025).

    # Arguments
    - `A::Array{precc,2}`                : 6x6 A matrix at radius r, which is used in the ODE for the solid-body problem.
    - `r::prec`                          : Radius at which to compute the A matrix.
    - `ρ::prec`                          : Density at radius r.
    - `g::prec`                          : Gravity at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `ω::prec`                          : Forcing frequency.
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
    function get_A!(A::Matrix, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        λ = Kd .- 2μ/3       # Lame's second param, which uses the drained compaction modulus
        S = ϕ/Kl + (α - ϕ)/K # Storavity, which uses liquid and solid grain bulk moduli  

        # First add the solid-body coefficients, but using drained moduli. 
        get_A!(A, r, ρ, g, μ, Kd, n; λ=λ)    # Note that here we replace the bulk modulus with the compaction modulus

        r_inv = 1.0/r
        β_inv = 1.0/(2μ + λ)

        # If there is a porous layer, now add the two-phase components
        if !iszero(ϕ)
            A[1,7] = α * β_inv

            A[3,1] += 1im * k*ρₗ^2 *g^2 * n*(n+1) / (ω*ηₗ) * r_inv^2
            A[3,5] += -(n+1)r_inv * 1im *(k*ρₗ^2*g*n)/(ω*ηₗ) * r_inv                               
            A[3,7] = 1im * (k*ρₗ*g*n*(n+1))/(ω*ηₗ)*r_inv^2 - 4μ*α*β_inv*r_inv 
            A[3,8] =  1im * (k*ρₗ^2*g^2*n*(n+1))/(ω*ηₗ)*r_inv^2 - 4ϕ*ρₗ*g*r_inv 
        
            A[4,7] = 2α*μ*r_inv * β_inv
            A[4,8] = ϕ*ρₗ*g*r_inv 
            
            A[5,8] = 4π*G*ρₗ*ϕ

            A[6,1] += -1im * 4π*G*n*(n+1)*r_inv * (k*ρₗ^2*g/(ω*ηₗ)*r_inv)
            A[6,5] = 1im*4π*n*(n+1)G*(ρₗ)^2*k*r_inv^2 / (ω*ηₗ)  
            A[6,7] = -1im *4π*n*(n+1)G*ρₗ*k*r_inv^2 / ( ω*ηₗ) 
            A[6,8] = 4π*G*(n+1)*r_inv * (ϕ*ρₗ - 1im * n*k*ρₗ^2*g/(ω*ηₗ)*r_inv) 
            
            A[7,1] = ρₗ*g*r_inv * ( 4 - 1im *(k*ρₗ*g*n*(n+1)/(ω*ϕ*ηₗ))*r_inv)  
            A[7,2] = -ρₗ*n*(n+1)*r_inv*g
            A[7,5] = -ρₗ*(n+1)r_inv * (1 - 1im*(k*ρₗ*g*n)/(ω*ϕ*ηₗ)*r_inv)  
            A[7,6] = ρₗ 
            A[7,7] = - 1im*(k*ρₗ*g*n*(n+1))/(ω*ϕ*ηₗ)*r_inv^2
            A[7,8] = -1im*ω*ϕ*ηₗ/k - 4π*G*(ρ - ϕ*ρₗ)*ρₗ + ρₗ*g*r_inv*(4 - 1im*(k*ρₗ*g*n*(n+1))/(ω*ϕ*ηₗ)*r_inv) 
        
            A[8,1] = r_inv*( 1im * k*ρₗ*g*n*(n+1)/(ω*ϕ*ηₗ)*r_inv - α/ϕ * 4μ*β_inv) 
            A[8,2] = α/ϕ * 2n*(n+1)*μ *β_inv * r_inv
            A[8,3] = -α/ϕ * β_inv 
            A[8,5] = -1im * k *ρₗ *n*(n+1) / (ω*ϕ*ηₗ)*r_inv^2 
            A[8,7] = 1im*k*n*(n+1)/(ω*ϕ*ηₗ)*r_inv^2 - 1/ϕ * (S + α^2 * β_inv) # If solid and liquid are compressible, keep the 1/M term
            A[8,8] = 1im * k *ρₗ*g *n*(n+1) / (ω*ϕ*ηₗ)*r_inv^2  - 2r_inv 
        end

    end


    """
        get_B(r1, r2, g1, g2, ρ, μ, K, n)

    Compute the 6x6 numerical integrator matrix, which integrates dy/dr from `r1` to `r2` for the solid-body problem.

    # Arguments
    - `r1::prec`                         : Starting radius for integration.
    - `r2::prec`                         : Ending radius for integration.
    - `g1::prec`                         : Gravity at radius r1.
    - `g2::prec`                         : Gravity at radius r2.
    - `ρ::prec`                          : Density at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `n::Int`                           : Tidal degree.

    # Returns
    - `B::Array{precc,2}`               : 6x6 numerical integrator matrix for integrating dy/dr from r1 to r2 for the solid-body problem.

    # Notes
    See 'get_B!' for definition.
    """ 
    function get_B(r1, r2, g1, g2, ρ, μ, K, n)
        B = zeros(precc, 6, 6)
        get_B!(B, r1, r2, g1, g2, ρ, μ, K, n)
        return B
    end


    """
        get_B!(B, r1, r2, g1, g2, ρ, μ, K)

    Compute the 6x6 numerical integrator matrix, which integrates dy/dr from `r1` to `r2` for the solid-body problem.
    `B` here represnts the RK4 integrator, given by Eq. S5.5 in Hay et al., (2025).

    # Arguments
    - `B::Array{precc,2}`                : 6x6 numerical integrator matrix for integrating dy/dr from r1 to r2 for the solid-body problem.
    - `r1::prec`                         : Starting radius for integration.
    - `r2::prec`                         : Ending radius for integration.
    - `g1::prec`                         : Gravity at radius r1.
    - `g2::prec`                         : Gravity at radius r2.
    - `ρ::prec`                          : Density at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `n::Int`                           : Tidal degree.

    # Notes
    See also [`get_B`](@ref)
    """
    function get_B!(B, r1, r2, g1, g2, ρ, μ, K, n)
        dr = r2 - r1
        rhalf = r1 + 0.5dr
        
        ghalf = g1 + 0.5*(g2 - g1)

        A1 = get_A(r1, ρ, g1, μ, K, n)
        Ahalf = get_A(rhalf, ρ, ghalf, μ, K, n)
        A2 = get_A(r2, ρ, g2, μ, K, n)
        
        k16 = zeros(precc, 6, 6)
        k26 = zeros(precc, 6, 6)
        k36 = zeros(precc, 6, 6)
        k46 = zeros(precc, 6, 6)

        k16 .= dr * A1 
        k26 .= dr * Ahalf * (I + 0.5k16)
        k36 .= dr * Ahalf * (I + 0.5k26)
        k46 .= dr * A2 * (I + k36) 

        # Only compute over the first six indices
        B[1:6,1:6] .= (I + 1.0/6.0 * (k16 + 2k26 + 2k36 + k46))
    end


    """
        get_B(r1, r2, g1, g2, ρ, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Compute the 8x8 numerical integrator matrix, which integrates dy/dr from `r1` to `r2` for the two-phase problem.

    # Arguments
    - `r1::prec`                         : Starting radius for integration.
    - `r2::prec`                         : Ending radius for integration.
    - `g1::prec`                         : Gravity at radius r1.
    - `g2::prec`                         : Gravity at radius r2.
    - `ρ::prec`                          : Density at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `ω::prec`                          : Forcing frequency.
    - `ρₗ::prec`                          : Liquid density at radius r.
    - `Kl::prec`                         : Liquid bulk modulus at radius r.
    - `Kd::prec`                         : Drained bulk modulus at radius r.
    - `α::prec`                          : Biot coefficient at radius r.
    - `ηₗ::prec`                          : Liquid viscosity at radius r.
    - `ϕ::prec`                          : Porosity at radius r.
    - `k::prec`                          : Permeability at radius r.
    - `n::Int`                           : Tidal degree.

    # Returns
    - `B::Array{precc,2}`               : 8x8 numerical integrator matrix for integrating dy/dr from r1 to r2 for the two-phase problem.

    # Notes
    See 'get_B!' for definition.
    """ 
    function get_B(r1, r2, g1, g2, ρ, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        B = zeros(precc, 8, 8)
        get_B!(B, r1, r2, g1, g2, ρ, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

        return B
    end


    """
        get_B!(B, r1, r2, g1, g2, ρ, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Compute the 8x8 numerical integrator matrix, which integrates dy/dr from `r1` to `r2` for the two-phase problem.
    `B` here represnts the RK4 integrator, given by Eq. S5.5 in Hay et al., (2025).

    # Arguments
    - `B::Array{precc,2}`                : 6x6 numerical integrator matrix for integrating dy/dr from r1 to r2 for the solid-body problem.
    - `r1::prec`                         : Starting radius for integration.
    - `r2::prec`                         : Ending radius for integration.
    - `g1::prec`                         : Gravity at radius r1.
    - `g2::prec`                         : Gravity at radius r2.
    - `ρ::prec`                          : Density at radius r.
    - `μ::prec`                          : Shear modulus at radius r.
    - `K::prec`                          : Bulk modulus at radius r.
    - `ω::prec`                          : Forcing frequency.
    - `ρₗ::prec`                          : Liquid density at radius r.
    - `Kl::prec`                         : Liquid bulk modulus at radius r.
    - `Kd::prec`                         : Drained bulk modulus at radius r.
    - `α::prec`                          : Biot coefficient at radius r.
    - `ηₗ::prec`                          : Liquid viscosity at radius r.
    - `ϕ::prec`                          : Porosity at radius r.
    - `k::prec`                          : Permeability at radius r.
    - `n::Int`                           : Tidal degree.

    # Notes
    See also [`get_B`](@ref)
    """
    function get_B!(B, r1, r2, g1, g2, ρ, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        dr = r2 - r1
        rhalf = r1 + 0.5dr
        
        ghalf = g1 + 0.5*(g2 - g1)

        Abot_p = zeros(precc, 8, 8)
        Amid_p = zeros(precc, 8, 8)
        Atop_p = zeros(precc, 8, 8)

        get_A!(Abot_p, r1, ρ, g1, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        get_A!(Amid_p, rhalf, ρ, ghalf, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        get_A!(Atop_p, r2, ρ, g2, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        
        k18 = zeros(precc, 8, 8)
        k28 = zeros(precc, 8, 8)
        k38 = zeros(precc, 8, 8)
        k48 = zeros(precc, 8, 8)

        k18 .= dr * Abot_p 
        k28 .= dr *  (Amid_p .+ 0.5Amid_p *k18) 
        k38 .= dr *  (Amid_p .+ 0.5*Amid_p *k28)
        k48 .= dr *  (Atop_p .+ Atop_p*k38) 

        I8 = SMatrix{8,8,precc}(I)

        B .= (I8 + 1.0/6.0 .* (k18 .+ 2*(k28 .+ k38) .+ k48))
    end


    """
        get_B_product!(Brod, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)

    Compute the product of the 8x8 B matrices within a primary layer. This is used to propgate the
    y solution across a single two-phase primary layer. Bprod is denoted by D in Eq. S5.14 in 
    Hay et al., (2025).

    # Arguments
    - `Bprod2::Array{precc,4}`           : 8x8x(nr-1)x(nlayers-1) array to store the B products across each secondary layer within each primary layer. 
    - `r::Array{prec,2}`                 : 2D array of layer boundaries.
    - `ρ::Array{prec,1}`                 : 1D array of layer densities. 
    - `g::Array{prec,2}`                 : 2D array of gravity values at the layer boundaries. 
    - `μ::Array{prec,1}`                 : 1D array of layer shear moduli.
    - `K::Array{prec,1}`                 : 1D array of layer bulk moduli.
    - `ω::prec`                          : Forcing frequency.
    - `ρₗ::Array{prec,1}`                 : 1D array of liquid densities at layer boundaries.
    - `Kl::Array{prec,1}`                : 1D array of liquid bulk moduli at layer boundaries.
    - `Kd::Array{prec,1}`                : 1D array of drained bulk moduli at layer boundaries.
    - `α::Array{prec,1}`                 : 1D array of Biot coefficients at layer boundaries.
    - `ηₗ::Array{prec,1}`                 : 1D array of liquid viscosities at layer boundaries.
    - `ϕ::Array{prec,1}`                 : 1D array of porosities at layer boundaries.
    - `k::Array{prec,1}`                 : 1D array of permeabilities at layer boundaries.
    - `n::Int`                           : Tidal degree.    
    """
    function get_B_product!(Bprod2, r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
        # Check dimensions of Bprod2

        nr = size(r)[1]

        Bstart = zeros(precc, 8, 8)
        B = zeros(precc, 8, 8)

        for i in 1:6
            Bstart[i,i] = 1
        end

        # if layer is porous, 
        # don't filter out y7 and y8 components
        if ϕ>0
            Bstart[7,7] = 1
            Bstart[8,8] = 1  
        end

        r1 = r[1]
        g1 = g[1]
        for j in 1:nr-1
            r2 = r[j+1]
            g2 = g[j+1]

            if ϕ>0 
                get_B!(B, r1, r2, g1, g2, ρ, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n)
            else
                get_B!(B, r1, r2, g1, g2, ρ, μ, K, n)
            end

            Bprod2[:,:,j] .= B * (j==1 ? Bstart : @view(Bprod2[:,:,j-1])) 

            r1 = r2
            g1 = g2 
        end
    end


    """
        compute_M(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; core="liquid")

    Compute the 4x4 M matrix, which relates the solution at the surface and porous layer interface to the core solution.
     
    # Arguments
    - `r::Array{prec,2}`                 : 2D array of layer boundaries.
    - `ρ::Array{prec,1}`                 : 1D array of layer densities. 
    - `g::Array{prec,2}`                 : 2D array of gravity values at the layer boundaries. 
    - `μ::Array{prec,1}`                 : 1D array of layer shear moduli.
    - `K::Array{prec,1}`                 : 1D array of layer bulk moduli.
    - `ω::prec`                          : Forcing frequency.
    - `ρₗ::Array{prec,1}`                 : 1D array of liquid densities at layer boundaries.
    - `Kl::Array{prec,1}`                : 1D array of liquid bulk moduli at layer boundaries.
    - `Kd::Array{prec,1}`                : 1D array of drained bulk moduli at layer boundaries.
    - `α::Array{prec,1}`                 : 1D array of Biot coefficients at layer boundaries.
    - `ηₗ::Array{prec,1}`                 : 1D array of liquid viscosities at layer boundaries.
    - `ϕ::Array{prec,1}`                 : 1D array of porosities at layer boundaries.
    - `k::Array{prec,1}`                 : 1D array of permeabilities at layer boundaries.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `core::String="liquid"`            : Type of core, either "liquid" or "solid". This is used to compute the starting vector for the numerical integration across the interior.

    # Returns
    - `M::Array{precc,2}`               : 4x4 M matrix, which is used to propagate the solution across the entire interior. 
    - `y1_4::Array{precc,4}`            : 4D array of the y solutions across each layer, which is used in the `compute_y` function to compute the solution vector across the interior.
    """
    function compute_M(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k, n; core="liquid")
        porous_layer = ϕ .> 0.0

        ## Convert parameters to the precision of precc:
        r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k = convert_params_to_prec(r, ρ, g, μ, K, ω, ρₗ, Kl, Kd, α, ηₗ, ϕ, k)

        sum(porous_layer) > 1.0 && error("Can only handle one porous layer for now!")

        nlayers = size(r)[2]
        nsublayers = size(r)[1]

        # Define starting vector as the core solution matrix, Y_r_C (Eq. S5.15)
        y_start = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], core, n; M=8, N=4)

        y1_4 = zeros(precc, 8,   4, nsublayers-1, nlayers)  # Four linearly independent y solutions
        
        for i in 2:nlayers
            Bprod = zeros(precc, 8, 8, nsublayers-1) # D matrix from Eq. S5.13
            @views get_B_product!(Bprod, r[:,i], ρ[i], g[:,i], μ[i], K[i], ω, ρₗ[i], Kl[i], Kd[i], α[i], ηₗ[i], ϕ[i], k[i], n)

            # Modify starting vector if the layer is porous
            # If a new porous layer (i.e., sitting on a non-porous layer)
            # reset the pore pressure and darcy flux
            if porous_layer[i] && !porous_layer[i-1]
                y_start[7,4] = 1.0          # Non-zero pore pressure
                y_start[8,4] = 0.0          # Zero radial Darcy flux
            elseif !porous_layer[i]
                y_start[7:8, :] .= 0.0      # Pore pressure and flux undefined
            end

            for j in 1:nsublayers-1
                y1_4[:,:,j,i] = @view(Bprod[:,:,j]) * y_start 
            end

            y_start[:,:] .= @view(y1_4[:,:,end,i])
        end

        # Get solution at the surface and porous layer interface (Eq. S5.20)
        M = zeros(precc, 4,4)
        
        M[1, :] .= y1_4[3,:,end,end] # Row 1 - Radial Stress
        M[2, :] .= y1_4[4,:,end,end] # Row 2 - Tangential Stress
        M[3, :] .= y1_4[6,:,end,end] # Row 3 - Potential Stress
        
        for i in 2:nlayers
            if porous_layer[i]
                M[4, :] .= y1_4[8,:,end,i]  #  Row 4 - Darcy flux (r = r_tp)
            end
        end
        
        return M, y1_4
    end


    """
        compute_y(r, g, M, R, y1_4, n; load=false)

    Compute the 8x1 solution vector at the surface and porous layer interface, which is used to compute the strain, 
    Darcy flux, and pore pressure at a particular radial level. This is given by Eq. S5.20 in Hay et al., (2025).

    # Arguments
    - `r::Array{prec,2}`                 : 2D array of layer boundaries.
    - `g::Array{prec,2}`                 : 2D array of gravity values at the layer boundaries. 
    - `M::Array{precc,2}`                : 4x4 M matrix, which is used to propagate the solution across the entire interior. 
    - `R::prec`                          : Surface radius of the body.
    - `y1_4::Array{precc,4}`             : 4D array of the y solutions across each layer, which is used in the `compute_y` function to compute the solution vector across the interior.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `load::Bool=false`                 : If true, compute the solution for a loaded problem.
    
    # Returns
    - `y::Array{ComplexF64,4}`           : 4D array of the solution vector y across the interior.
    """
    function compute_y(r, g, M, R, y1_4, n; load=false)

        nlayers = size(r)[2]
        nsublayers = size(r)[1]

        b = zeros(precc, 4)
        if load
            b[1] = -(2n+1)*g[end,end]/(4π*(R)^2)
            b[3] = -(2n+1)*G/(R)^2
        else
            b[3] = (2n+1)/R 
        end

        C = M \ b

        y = zeros(ComplexF64, 8, nsublayers-1, nlayers)

        for i in 2:nlayers
            for j in 1:nsublayers-1
                y[:,j,i] = @view(y1_4[:,:,j,i])*C
            end
        end

        return y
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

        @views Y    = solid1d_mush.Y[i,:,:]
        @views dYdθ = solid1d_mush.dYdθ[i,:,:]
        @views dYdϕ = solid1d_mush.dYdϕ[i,:,:]
        @views Z    = solid1d_mush.Z[i,:,:]
        @views X    = solid1d_mush.X[i,:,:]

        y1 = y[1]
        y2 = y[2]
        y3 = y[3]
        y4 = y[4]
        y7 = y[7]

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

        @views Y    = solid1d_mush.Y[i,:,:]
        @views dYdθ = solid1d_mush.dYdθ[i,:,:]
        @views dYdϕ = solid1d_mush.dYdϕ[i,:,:]
        
        y1 = y[1]
        y5 = y[5]
        y7 = y[7]
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

        @views Y    = solid1d_mush.Y[i,:,:]
        @views dYdθ = solid1d_mush.dYdθ[i,:,:]
        @views dYdϕ = solid1d_mush.dYdϕ[i,:,:]
        
        y7 = y[7]
        
        p[:,:] = Y * y7 
    end


    """
        get_heating_profile(y, r, ρ, g, μ, Ks, ω, ρl, Kl, Kd, α, ηl, ϕ, k, n; lay=nothing)

    Get the radial volumetric heating for two-phase tides and eccentricity forcing,
    assuming synchronous rotation.
    Heating rate is computed with numerical integration using the solution 
    `y` returned by [`compute_y`](@ref), 
    using Eq. 2.39a/b/c integrated over solid angle. 
    The heating profile for a specific layer is specified with `lay`, otherwise all 
    layers will be caclulated.

    # Arguments
    - `y::Array{ComplexF64,4}`           : 4D array of the solution vector y across the interior, returned by `compute_y`.
    - `r::Array{Float64,2}`              : 2D array of layer boundaries.
    - `ρ::Array{Float64,1}`              : 1D array of layer densities.
    - `g::Array{Float64,2}`              : 2D array of gravity values at the layer boundaries.
    - `μ::Array{Float64,1}`              : 1D array of layer shear moduli.
    - `Ks::Array{Float64,1}`              : 1D array of layer bulk moduli for shear dissipation.
    - `ω::Float64`                       : Tidal frequency in radians per second.
    - `ρl::Array{Float64,1}`             : 1D array of liquid densities at layer boundaries.
    - `Kl::Array{Float64,1}`             : 1D array of liquid bulk moduli at layer boundaries.
    - `Kd::Array{Float64,1}`             : 1D array of drained bulk moduli at layer boundaries.
    - `α::Array{Float64,1}`              : 1D array of Biot coefficients at layer boundaries.
    - `ηl::Array{Float64,1}`             : 1D array of liquid viscosities at layer boundaries.
    - `ϕ::Array{Float64,1}`             : 1D array of porosities at layer boundaries.
    - `k::Array{Float64,1}`              : 1D array of permeabilities at layer boundaries.
    - `n::Int`                           : Tidal degree.

    # Keyword Arguments
    - `lay::Int=nothing`                 : If specified, compute the heating profile for only the layer corresponding to this index. Otherwise, compute for all layers.

    # Returns
    - `Eμ_tot::Array{Float64,1}`         : 1D array of total power dissipated in each primary layer due to shear, in W.
    - `Eμ_vol::Array{Float64,2}`         : 2D array of angular averaged volumetric heating profiles in W/m^3 for dissipation due to shear, with dimensions corresponding to the secondary layer and primary layer indices.
    - `Eκ_tot::Array{Float64,1}`         : 1D array of total power dissipated in each primary layer due to compaction, in W.
    - `Eκ_vol::Array{Float64,2}`         : 2D array of angular averaged volumetric heating profiles in W/m^3 for dissipation due to compaction, with dimensions corresponding to the secondary layer and primary layer indices.
    - `El_tot::Array{Float64,1}`         : 1D array of total power dissipated in each primary layer due to Darcy flow, in W.
    - `El_vol::Array{Float64,2}`         : 2D array of angular averaged volumetric heating profiles in W/m^3 for dissipation due to Darcy flow, with dimensions corresponding to the secondary layer and primary layer indices.
    """
    function get_heating_profile(y, r, ρ, g, μ, Ks, ω, ρl, Kl, Kd, α, ηl, ϕ, k, n; lay=nothing)
        dres = deg2rad(solid1d_mush.res)

        @views clats = solid1d_mush.clats[:]
        @views lons = solid1d_mush.lons[:]

        nsublay = size(r)[1]
        nlay = size(r)[2]
        nlats = length(clats)
        nlons = length(lons)

        ϵ = zeros(ComplexF64, nlats, nlons, 6, nsublay-1, nlay)
        d_disp = zeros(ComplexF64, nlats, nlons, 3, nsublay-1, nlay)
        p = zeros(ComplexF64, nlats, nlons, nsublay-1, nlay)

        ϵs = zero(ϵ)
        d_disps = zero(d_disp)
        ps = zero(p)

        # # Retrieve stain tensor 
        x = 1
        @views Y    = solid1d_mush.Y[x,:,:]

        for i in 2:nlay # Loop of layers
            ρr = ρ[i]
            Ksr = Ks[i]
            μr = μ[i]
            ρlr = ρl[i]
            Klr = Kl[i]
            Kdr = Kd[i]
            αr = α[i]
            ηlr = ηl[i]
            ϕr = ϕ[i]
            kr = k[i]

            for j in 1:nsublay-1 # Loop over sublayers 
                @views yrr = y[:,j,i]
                (y1, y2, y3, y4, y5, y6, y7, y8) = yrr
                
                rr = r[j,i]
                gr = g[j,i]

                compute_strain_ten!(@view(ϵ[:,:,:,j,i]), yrr, n, rr, ρr, gr, μr, Ksr, ω, ρlr, Klr, Kdr, αr, ηlr, ϕr, kr)
                
                if ϕ[i] > 0 
                    compute_darcy_displacement!(@view(d_disp[:,:,:,j,i]), yrr, n, rr, ω, ϕr, ηlr, kr, gr, ρlr)
                    compute_pore_pressure!(@view(p[:,:,j,i]), yrr, n)
                end

                p[:,:,j,i] .= y7 * Y    # pore pressure

            end
        end

        ϵs .+= ϵ
        d_disps .+= d_disp
        ps .+= p

        # Shear heating in the solid
        Eμ = zeros(  (size(ϵ)[1], size(ϵ)[2], size(ϵ)[4], size(ϵ)[5]) )
        Eμ_tot = zeros(  (nlay) )
        Eμ_vol = zeros(  size(r) )

        # Darcy dissipation in the liquid
        El = zeros(  (size(ϵ)[1], size(ϵ)[2], size(ϵ)[4], size(ϵ)[5]) )
        El_tot = zeros(  (nlay) )
        El_vol = zeros(  size(r) )

        # Bulk dissipation in the solid
        Eκ = zero(Eμ)
        Eκ_tot = zero( Eμ_tot )
        Eκ_vol = zero( Eμ_vol )


        if isnothing(lay)
            rstart = 2
            rend = nlay
        else
            rstart = lay
            rend = lay
        end

        for j in rstart:rend   # loop from CMB to surface          
            for i in 1:nsublay-1
                dr = (r[i+1, j] - r[i, j])
                dvol = 4π/3 * (r[i+1, j]^3 - r[i, j]^3)

                @views Eμ[:,:,i, j] = sum(abs.(ϵs[:,:,1:3,i,j]).^2, dims=3) .+ 
                                    2sum(abs.(ϵs[:,:,4:6,i,j]).^2, dims=3) .- 
                                    1/3 .* abs.(sum(ϵs[:,:,1:3,i,j], dims=3)).^2
                Eμ[:,:,i, j] .*= ω * imag(μ[j])

                @views Eκ[:,:,i, j] = ω/2 *imag(Kd[j]) * abs.(sum(ϵs[:,:,1:3,i,j], dims=3)).^2
                if ϕ[j] > 0
                    @views Eκ[:,:,i, j] .+= ω/2 *imag(Kd[j]) * (abs.(ps[:,:,i,j]) ./ Ks[j]).^2
                end

                # Integrate across r to find dissipated energy per unit area
                @views Eμ_vol[i,j] = sum(sin.(clats) .* (Eμ[:,:,i,j])  * dres^2) * r[i,j]^2.0 * dr / dvol

                @views Eκ_vol[i,j] = sum(sin.(clats) .* (Eκ[:,:,i,j])  * dres^2) * r[i,j]^2.0 * dr / dvol
        
                if ϕ[j] > 0            
                    @views El[:,:,i, j] = 0.5 *  ϕ[j]^2 * ω^2 * ηl[j]/k[j] * (abs.(d_disps[:,:,1,i,j]).^2 + abs.(d_disps[:,:,2,i,j]).^2 + abs.(d_disps[:,:,3,i,j]).^2)
                    @views El_vol[i,j] = sum(sin.(clats) .* (El[:,:,i,j])  * dres^2) * r[i,j]^2.0 * dr / dvol

                end

            end
            
            # Average over sublayers to get layer volumetric heating (W/m³)
            Eμ_tot[j] = sum(Eμ_vol[1:nsublay-1,j] .* diff(r[1:nsublay,j])) / (r[nsublay,j] - r[1,j])
            Eκ_tot[j] = sum(Eκ_vol[1:nsublay-1,j] .* diff(r[1:nsublay,j])) / (r[nsublay,j] - r[1,j])
            El_tot[j] = sum(El_vol[1:nsublay-1,j] .* diff(r[1:nsublay,j])) / (r[nsublay,j] - r[1,j])

        end

        return (Eμ_tot, Eμ_vol), (Eκ_tot, Eκ_vol), (El_tot, El_vol)
    end

end
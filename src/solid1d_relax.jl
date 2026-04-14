

module solid1d_relax
    
    using LinearAlgebra
    using DoubleFloats
    using AssociatedLegendrePolynomials    
    using StaticArrays
    using SpecialFunctions
    using SparseArrays
    using Statistics

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
        convert_params_to_prec(r, ρ, g, μ, κs)

    Convert input parameters into the required precision.
    # Arguments
    - `r::Array{Float64,2}`               : 2D array of layer boundaries.
    - `ρ::Array{Float64,1}`               : 1D array of layer densities. 
    - `g::Array{Float64,2}`               : 2D array of gravity values at the layer boundaries. 
    - `μ::Array{Float64,1}`               : 1D array of layer shear moduli.
    - `κs::Array{Float64,1}`              : 1D array of layer bulk moduli. 
    
    # Returns
    Tuple of converted parameters in the required precision.
    """
    function convert_params_to_prec(r, ρ, g, μ, κs)
        r_prec = convert(Array{prec}, r)
        ρ_prec = convert(Array{prec}, ρ)
        g_prec = convert(Array{prec}, g)
        μ_prec = convert(Array{precc}, μ)
        κs_prec = convert(Array{precc}, κs)

        return (r_prec,  ρ_prec, g_prec, μ_prec, κs_prec)
    end


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


    function safe_sph_besselj(l, x)
        try
            return sphericalbesselj(l, x)
        catch
            # asymptotic fallback
            return sin(x - l*pi/2) / x
        end
    end


    """
        get_Ic(ω, r, ρ, g, μ, K, type, n; M=6, N=3)
            
    Get the core solution vector.
    
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

            jl_n   = safe_sph_besselj(n, x64)
            jl_np1 = safe_sph_besselj(n+1, x64)

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
        else # incompressible solid core
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
    function get_A(ω, r, ρ, g, μ, K, n; G0=1, λ=nothing)
        A = zeros(precc, 6, 6) 
        get_A!(A, ω, r, ρ, g, μ, K, n; G0=G0, λ=λ)
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


    function solve_radial_system(r, ρ, g, μ, K, ω, n, R_planet; core="inertial")

        Nr = length(r)

        # scaling
        # R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(1, 1, 1)
        R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(9.8e6, 4.3e22, 4.3e3)

        # Extract sub-blocks for BC scaling
        Su = S[1:3, 1:3]
        Sl = S[4:6, 4:6]
        Sinv_u = Sinv[1:3, 1:3]
        Sinv_l = Sinv[4:6, 4:6]

        ωs = ω * s0
        rs = r ./ R0
        ρs = ρ ./ ρ0
        gs = g ./ g0
        μs = μ ./ μ0
        Ks = K ./ μ0

        # boundary conditions
        B1 = zeros(precc, 3, 6)

        function isbad(A)
            any(x -> !isfinite(real(x)) || !isfinite(imag(x)), A)
        end

        try
            B1 .= get_core_bc!(ω, r[1], ρ[1], g[1], μ[1], K[1], core, n; G0=1)

            if isbad(B1)
                println("NaN/Inf detected in core BC → switching to liquid core")
                B1 .= get_core_bc!(ω, r[1], ρ[1], g[1], μ[1], K[1], "liquid", n; G0=1)
            end

        catch e
            println("Error in core BC → switching to liquid core")
            B1 .= get_core_bc!(ω, r[1], ρ[1], g[1], μ[1], K[1], "liquid", n; G0=1)
        end
        
        BN, b = get_surface_bc!(R_planet, n)
        # BN, b = get_surface_bc!(R_planet, n)

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

    # function solve_radial_system(r, ρ, g, μ, K, ω, n, R_planet; core="inertial")

    #     Nr = length(r)

    #     # scaling
    #     # R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(1, 1, 1)
    #     R0, M0, s0, ρ0, G0, g0, μ0, S, Sinv = get_scales(9.8e6, 4.3e22, 4.3e3)

    #     ωs = ω * s0
    #     rs = r ./ R0
    #     ρs = ρ ./ ρ0
    #     gs = g ./ g0
    #     μs = μ ./ μ0
    #     Ks = K ./ μ0

    #     println(g0)

    #     # boundary conditions
    #     B1 = zeros(precc, 3, 6)

    #     # function isbad(A)
    #     #     any(x -> !isfinite(real(x)) || !isfinite(imag(x)), A)
    #     # end

    #     # try
    #     #     B1 .= get_core_bc!(ωs, rs[1], ρs[1], gs[1], μs[1], Ks[1], core, n)

    #     #     if isbad(B1)
    #     #         println("NaN/Inf detected in core BC → switching to liquid core")
    #     #         B1 .= get_core_bc!(ωs, rs[1], ρs[1], gs[1], μs[1], Ks[1], "liquid", n)
    #     #     end

    #     # catch e
    #     #     println("Error in core BC → switching to liquid core")
    #     #     B1 .= get_core_bc!(ωs, rs[1], ρs[1], gs[1], μs[1], Ks[1], "liquid", n)
    #     # end
    #     B1 .= get_core_bc!(ωs, rs[1], ρs[1], gs[1], μs[1], Ks[1], core, n; G0=G0)
    #     BN, b = get_surface_bc!(R_planet/R0, n)

    #     # storage
    #     R = Vector{Matrix{precc}}(undef, Nr)

    #     # first layer (n = 1)
    #     dr = rs[2] - rs[1]

    #     A1 = get_A(ωs, rs[1], ρs[1], gs[1], μs[1], Ks[1], n; G0=G0)
    #     A2 = get_A(ωs, rs[2], ρs[2], gs[2], μs[2], Ks[2], n; G0=G0)

    #     I6 = Matrix{precc}(I, 6, 6)

    #     C1 =  I6 + 0.5 * dr * A1
    #     D2 = -I6 + 0.5 * dr * A2

    #     # split matrices
    #     C1u, C1l = C1[1:3, :], C1[4:6, :]
    #     D2u, D2l = D2[1:3, :], D2[4:6, :]

    #     # build S1 and Q1
    #     S1 = [B1; C1u]              # 6×6
    #     Q1 = [zeros(3,6); D2u]      # 6×6

    #     # initial recursion
    #     R[1] = -pinv(S1) * Q1

    #     # forward recursion
    #     for i in 2:Nr-1

    #         dr = rs[i+1] - rs[i]

    #         A_n  = get_A(ωs, rs[i],   ρs[i],   gs[i],   μs[i],   Ks[i],   n; G0=G0)
    #         A_np = get_A(ωs, rs[i+1], ρs[i+1], gs[i+1], μs[i+1], Ks[i+1], n; G0=G0)

    #         Cn  =  I6 + 0.5 * dr * A_n
    #         Dnp = -I6 + 0.5 * dr * A_np

    #         # split
    #         Cn_u, Cn_l = Cn[1:3, :], Cn[4:6, :]
    #         Dnp_u, Dnp_l = Dnp[1:3, :], Dnp[4:6, :]

    #         # build blocks
    #         Pn = [Cn_l; zeros(3,6)]
    #         Sn = [Dnp_l; Cn_u]
    #         Qn = [zeros(3,6); Dnp_u]

    #         # recursion
    #         Xn = Pn * R[i-1] + Sn
    #         R[i] = -pinv(Xn) * Qn
    #     end

    #     # final layer (n = N)
    #     dr = rs[Nr] - rs[Nr-1]


    #     A_Nm = get_A(ωs, rs[end-1], ρs[end-1], gs[end-1], μs[end-1], Ks[end-1], n; G0=G0)
    #     A_N  = get_A(ωs, rs[end],   ρs[end],   gs[end],   μs[end],   Ks[end],   n; G0=G0)

    #     CNm =  I6 + 0.5 * dr * A_Nm
    #     DN  = -I6 + 0.5 * dr * A_N

    #     CNm_l = CNm[4:6, :]
    #     DN_l  = DN[4:6, :]

    #     PN = [CNm_l; zeros(3,6)]
    #     SN = [DN_l; BN]

    #     XN = PN * R[Nr-1] + SN

    #     return XN, R, b
    # end


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


    function get_surface_bc!(R, n)
        
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


    function compute_y_relaxation(r, ρ, g, μ, K, ω, n, R)

        if ω < 1e-6
            core = "liquid"
        else
            core = "inertial"
        end

        # core = "liquid"

        XN, R, b = solve_radial_system(r, ρ, g, μ, K, ω, n, R; core=core)

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
    - `Eμ_vol::Array{Float64,2}`         : 2D array of angular averaged volumetric heating profiles in W/m^3 for dissipation due to shear, with dimensions corresponding to the secondary layer and primary layer indices.
    - `Eκ_tot::Array{Float64,1}`         : 1D array of total power dissipated in each primary layer due to compaction, in W.
    - `Eκ_vol::Array{Float64,2}`         : 2D array of angular averaged volumetric heating profiles in W/m^3 for dissipation due to compaction, with dimensions corresponding to the secondary layer and primary layer indices.
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

            Eμ_vol[i] = sum(weight .* Eμ_loc * dres^2) * rr^2 * dr / dvol
            Eκ_vol[i] = sum(weight .* Eκ_loc * dres^2) * rr^2 * dr / dvol

            # accumulate totals
            Eμ_tot[i] = Eμ_vol[i] 
            Eκ_tot[i] = Eκ_vol[i] 
        end

        return (Eμ_tot, Eμ_vol), (Eκ_tot, Eκ_vol)
    end

end
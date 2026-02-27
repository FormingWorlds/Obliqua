

module solid1d
    
    using LinearAlgebra
    using DoubleFloats
    using AssociatedLegendrePolynomials    
    using StaticArrays


    prec = BigFloat
    precc = Complex{BigFloat}

    const G::prec       = prec(6.6743e-11)       # m^3 kg^-1 s^-2

    M = 6

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


    function Ynm(n, m, theta, phi)
        return Plm.(n, m, cos.(theta)) .* exp.(1im * m .* phi)
    end


    """
        define_spherical_grid(res)

    Create the spherical grid of angular resolution `res` in degrees. This is used for 
    numerical integrations over solid angle. A new grid can easily be defined by 
    recalling the function with a new `res`.

    The grid is internal to solid1d, but can be accessed with 

        solid1d.clats[:] # colatitude grid
        solid1d.lons[:]  # longitude grid
    """
    function define_spherical_grid(res, n, m)
        solid1d.res = res

        # θ and φ grids
        lons = deg2rad.(collect(0:res:360-0.001))'
        clats = deg2rad.(collect(0:res:180))
        clats[1] += 1e-6
        clats[end] -= 1e-6

        # allocate arrays
        solid1d.Y    = zeros(ComplexF64, 1, length(clats), length(lons))
        solid1d.dYdθ = similar(solid1d.Y)
        solid1d.dYdϕ = similar(solid1d.Y)
        solid1d.Z    = similar(solid1d.Y)
        solid1d.X    = similar(solid1d.Y)

        sinθ = sin.(clats)
        cosθ = cos.(clats)
        cotθ = cosθ ./ sinθ
        cscθ = csc.(clats)

        # Normalization factor for spherical harmonics
        norm = sqrt((2*n+1)  * factorial(n-m) / (4π * factorial(n+m)))
        
        i = 1

        # Y
        solid1d.Y[i,:,:] = Ynm(n,m,clats,lons)

        # ∂Y/∂θ
        Pn = Plm.(n, m, cosθ)
        if n > m
            Pn_1 = Plm.(n-1, m, cosθ)
            dPdθ = (n .* cosθ .* Pn .- (n + m) .* Pn_1) ./ (sinθ)
        else
            # m == n -> P_{n-1}^m = 0
            dPdθ = (n .* cosθ .* Pn) ./ (sinθ)
        end
        solid1d.dYdθ[i,:,:] .= dPdθ .* exp.(1im .* m .* lons)

        # ∂Y/∂ϕ
        solid1d.dYdϕ[i,:,:] .= 1im * m .* solid1d.Y[i,:,:]

        # Z = 2 ((1/sinθ) ∂²Y/∂θ∂ϕ - cotθ cscθ ∂Y/∂ϕ)
        solid1d.Z[i,:,:] .= 2 .* (1im * m ./ sinθ .* solid1d.dYdθ[i,:,:] .- cotθ .* cscθ .* solid1d.dYdϕ[i,:,:])

        # X = -2 (cotθ ∂Y/∂θ + csc²θ ∂²Y/∂ϕ²) - n(n+1)) Y
        solid1d.X[i,:,:] .= -2 .* (cotθ .* solid1d.dYdθ[i,:,:] .- cscθ.^2 .* m^2 .* solid1d.Y[i,:,:]) .- n*(n+1) .* solid1d.Y[i,:,:]

        # Normalize
        solid1d.Y[i,:,:]    .*= norm
        solid1d.dYdθ[i,:,:] .*= norm
        solid1d.dYdϕ[i,:,:] .*= norm
        solid1d.Z[i,:,:]    .*= norm
        solid1d.X[i,:,:]    .*= norm

        # save grids
        solid1d.clats = clats
        solid1d.lons  = lons
    end


    "Convert input parameters into the required precision."
    function convert_params_to_prec(r, ρ, g, μ, κs)
        r_prec = convert(Array{prec}, r)
        ρ_prec = convert(Array{prec}, ρ)
        g_prec = convert(Array{prec}, g)
        μ_prec = convert(Array{precc}, μ)
        κs_prec = convert(Array{precc}, κs)

        return (r_prec,  ρ_prec, g_prec, μ_prec, κs_prec)
    end


    "Get the core solution vector"
    function get_Ic(r, ρ, g, μ, type, n; M=6, N=3)
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

    See also [`get_A!`](@ref)
    """
    function get_A(r, ρ, g, μ, K, n)
        A = zeros(precc, 6, 6) 
        get_A!(A, r, ρ, g, μ, K, n)
        return A
    end


    """
        get_A!(r, ρ, g, μ, K, n; λ=nothing)

    Compute the 6x6 `A` matrix in the ODE for the solid-body problem. These correspond to 
    the coefficients given in Equation S4.6 in Hay et al., (2025) when α=φ=0, as well as Sabadini and Vermeersen 
    (2016) Eq. 1.95.

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


    "See 'get_B!' for definition."
    function get_B(r1, r2, g1, g2, ρ, μ, K, n)
        B = zeros(precc, 6, 6)
        get_B!(B, r1, r2, g1, g2, ρ, μ, K, n)
        return B
    end


    """
        get_B!(B, r1, r2, g1, g2, ρ, μ, K)

    Compute the 6x6 numerical integrator matrix, which integrates dy/dr from `r1` to `r2` for the solid-body problem.

    `B` here represnts the RK4 integrator, given by Eq. S5.5 in Hay et al., (2025).

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

        I6 = SMatrix{6,6,precc}(I)

        # Only compute over the first six indices
        B[1:6,1:6] .= (I6 + 1.0/6.0 * (k16 + 2k26 + 2k36 + k46))
    end


    """
        get_B_product!(Brod, r, ρ, g, μ, K)

    Compute the product of the 6x6 B matrices within a primary layer. This is used to propgate the
    y solution across one single-phase (solid) primary layer.

    Bprod is denoted by D in Eq. S5.14 in Hay et al., (2025).
    """
    function get_B_product!(Bprod2, r, ρ, g, μ, K, n)
        Bstart = zeros(precc, 6, 6)
        B = zeros(precc, 6, 6)

        for i in 1:6
            Bstart[i,i,1] = 1.0
        end

        nr = size(r)[1]

        r1 = r[1]
        for j in 1:nr-1
            r2 = r[j+1]
            g1 = g[j]
            g2 = g[j+1]

            get_B!(B, r1, r2, g1, g2, ρ, μ, K, n)
            Bprod2[:,:,j] .= B * (j==1 ? Bstart : Bprod2[:,:,j-1])

            r1 = r2
        end
    end


    function compute_M(r, ρ, g, μ, K, n; core="liquid", load=false)
        r, ρ, g, μ, K = convert_params_to_prec(r, ρ, g, μ, K)

        nlayers = size(r)[2]
        nsublayers = size(r)[1]

        y_start = get_Ic(r[end,1], ρ[1], g[end,1], μ[1], core, n; M=6, N=3)

        y1_4 = zeros(precc, 6, 3, nsublayers-1, nlayers) # Three linearly independent y solutions
                
        for i in 2:nlayers
            Bprod = zeros(precc, 6, 6, nsublayers-1)
            @views get_B_product!(Bprod, r[:, i], ρ[i], g[:, i], μ[i], K[i], n)

            for j in 1:nsublayers-1
                y1_4[:,:,j,i] = @view(Bprod[:,:,j]) * y_start 
            end

            y_start[:,:] .= @view(y1_4[:,:,end,i])   # Set starting vector for next layer
        end

        M = zeros(precc, 3,3)

        M[1, :] .= y1_4[3,:,end,end]  # Row 1 - Radial Stress 
        M[2, :] .= y1_4[4,:,end,end]  # Row 2 - Tangential Stress
        M[3, :] .= y1_4[6,:,end,end]  # Row 3 - Potential Stress
        
        return M, y1_4
    end


    function compute_y(r, g, M, R, y1_4, n; load=false)

        nlayers = size(r)[2]
        nsublayers = size(r)[1]

        b = zeros(precc, 3)
        if load
            b[1] = -(2n+1)*g[end,end]/(4π*(R)^2) ## Fix g[end,end] !!!!!
            b[3] = -(2n+1)*G/(R)^2
        else
            b[3] = (2n+1)/R 
        end

        C = M \ b

        y = zeros(ComplexF64, 6, nsublayers-1, nlayers)

        for i in 2:nlayers
            for j in 1:nsublayers-1
                y[:,j,i] = @view(y1_4[:,:,j,i])*C
            end
        end

        return y
    end


    """
        compute_strain_ten!(ϵ, y, m, rr, ρr, gr, μr, Ksr, ω, ρlr, Klr, Kdr, αr, ηlr, ϕr, kr)

        Calculate the strain tensor ϵ at a particular radial level. 
    """
    function compute_strain_ten!(ϵ, y, n, rr, ρr, gr, μr, Ksr)
        i = 1

        @views Y    = solid1d.Y[i,:,:]
        @views dYdθ = solid1d.dYdθ[i,:,:]
        @views dYdϕ = solid1d.dYdϕ[i,:,:]
        @views Z    = solid1d.Z[i,:,:]
        @views X    = solid1d.X[i,:,:]

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
        function get_heating_profile(y, r, ρ, g, μ, κ, ω, ecc; lay=nothing)

    Get the radial volumetric heating for solid-body tides and eccentricity forcing,
    assuming synchronous rotation.
    Heating rate is computed with numerical integration using the solution 
    `y` returned by [`compute_y`](@ref), 
    using Eq. 2.39a/b integrated over solid angle. 
    The forcing magnitude is determined by orbital frequency `ω` and eccentricity `ecc`. 
    The heating profile for a specific layer is specified with `lay`, otherwise all 
    layers will be caclulated.

    Returns the angular averaged volumetric heating profiles in W/m^3 for dissipation due to shear 
    and compaction, as well as the power dissipated in each primary layer in W/m^3.
    """
    function get_heating_profile(y, r, ρ, g, μ, κ, n, ω; lay=nothing)

        dres = deg2rad(solid1d.res)

        @views clats = solid1d.clats[:]
        @views lons  = solid1d.lons[:]

        nsublay = size(r)[1]
        nlay = size(r)[2]
        nlats = length(clats)
        nlons = length(lons)

        ϵ       = zeros(ComplexF64, nlats, nlons, 6, nsublay, nlay)
        ϵs      = zero(ϵ)

        # Retrieve stain tensor 
        for i in 2:nlay # Loop of layers
            ρr = ρ[i]
            κr = κ[i]
            μr = μ[i]

            for j in 1:nsublay-1 # Loop over sublayers 
                @views yrr = y[1:6,j,i]
                
                rr = r[j,i]
                gr = g[j,i]

                compute_strain_ten!(@view(ϵ[:,:,:,j,i]), yrr, n, rr, ρr, gr, μr, κr)
            end
        end

        ϵs .+= ϵ 

        Eμ      = zeros(  (nlats, nlons, nsublay-1, nlay) )

        Eμ_tot = zeros(  (nlay) )     # Total heating rate (W)
        Eμ_vol = zeros(  size(r) )    # Spherical avg. volumetric heating rate in each layer (W/m^3)

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

        # Compute dissipation (volumetric heating in W/m^3)
        for j in rstart:rend    
            for i in 1:nsublay-1
                dr = r[i+1, j] - r[i, j]                  # radial thickness of sublayer
                dvol = 4π/3 * (r[i+1, j]^3 - r[i, j]^3)   # volume of spherical shell

                Eμ[:,:,i, j] = sum(abs.(ϵs[:,:,1:3,i,j]).^2, dims=3) .+
                                2sum(abs.(ϵs[:,:,4:6,i,j]).^2, dims=3) .-
                                1/3 .* abs.(sum(ϵs[:,:,1:3,i,j], dims=3)).^2
                Eμ[:,:,i, j] .*= ω * imag(μ[j])

                Eκ[:,:,i, j] = ω/2 * imag(κ[j]) * abs.(sum(ϵs[:,:,1:3,i,j], dims=3)).^2

                # Angular integration to get volumetric heating
                Eμ_vol[i,j] = sum(sin.(clats) .* Eμ[:,:,i,j] * dres^2) * r[i,j]^2.0 * dr / dvol
                Eκ_vol[i,j] = sum(sin.(clats) .* Eκ[:,:,i,j] * dres^2) * r[i,j]^2.0 * dr / dvol
                
            end

            # Average over sublayers to get layer volumetric heating (W/m^3)
            Eμ_tot[j] = sum(Eμ_vol[1:nsublay-1,j] .* diff(r[1:nsublay,j])) / (r[nsublay,j] - r[1,j])
            Eκ_tot[j] = sum(Eκ_vol[1:nsublay-1,j] .* diff(r[1:nsublay,j])) / (r[nsublay,j] - r[1,j])
        end


        return (Eμ_tot, Eμ_vol), (Eκ_tot, Eκ_vol) 
    end

end
# Helper functions shared by the Obliqua models: argument checks, the tidal work of the energy closure, and diagnostics
# of a solved column (face residuals, boundary-layer resolution). Included as its own module (`using .helper`) by the
# modules that need it, like constants.jl and common.jl.

module helper

    include("constants.jl")
    using .constants

    export check_arg, check_option, tidal_work, condition_label, face_residuals, boundary_layer_resolution


    """
        check_arg(cond, msg)

    Throw `ArgumentError(msg)` unless `cond` holds. For checks on the arguments of a public function, e.g.
    `check_arg(dr_max > dr_min > 0, "column_spec: need dr_max > dr_min > 0, got \$dr_min, \$dr_max")`.
    """
    check_arg(cond::Bool, msg::AbstractString) = cond || throw(ArgumentError(msg))

    """
        check_option(fn, name, value, options)

    Throw an `ArgumentError` naming the function `fn`, the argument `name` and the allowed `options` unless `value` is
    one of them, e.g. `check_option("column_spec", "slip", slip, (:freeslip, :noslip))`.
    """
    check_option(fn::AbstractString, name::AbstractString, value, options) =
        value in options || throw(ArgumentError("$fn: $name must be one of $(join(repr.(options), ", ")), got $(repr(value))"))


    """
        tidal_work(k, n, r, ω)

    Tidal work per |U_T|^2 [W] for a Love number k of degree n referenced at radius r [m]:
    (2n+1) r w / (8 pi G) (-Im k). Equals the integrated dissipation for an energy-consistent solution.
    """
    tidal_work(k, n, r, ω) = (2n + 1) * r * ω / (8π * G) * (-imag(k))


    """
        boundary_layer_resolution(cs)

    Poroelastic boundary-layer thickness of every mush node over its cell size (appendix eq:delta):

        delta = sqrt(2 k |M_B| / (eta_l omega)),   1/M_B = S + alpha^2 / beta_d,   beta_d = K_d + 4 mu / 3,
        S = phi / K_l + (alpha - phi) / K_s.

    The melt flux W vanishes at an impermeable face and recovers over delta, and at a porosity step the pore pressure
    and W relax over delta on each side. A collapsed mush (K_d = 0, mu ~ G_inf) has |M_B| ~ beta_d and a very thin
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
        condition_label(c, lo, hi)

    Readable form of a condition, e.g. "U⁻ + W⁻ = U⁺", "Y⁻ = 0".
    """
    function condition_label(c, lo::Symbol, hi::Symbol)
        name(s) = ("U", "V", "Φ", "P", "X", "Y", "Ψ", "W")[s]     # solver slots: slot 8 is the melt flux W
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

end

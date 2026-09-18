# src/fourier_tools.jl

# ==============================================================================
# 3D SPHERICAL FOURIER TRANSFORMS AND S(q) <-> g(r) INVERSIONS
# ==============================================================================
# This module implements continuous 3D isotropic Fourier sine transformations
# allowing conversion between the real-space pair distribution function g(r)
# and the reciprocal-space static structure factor S(q).
#
# References:
# - Hansen, J.-P., & McDonald, I. R. (2013). Theory of Simple Liquids (4th ed.). Academic Press.
# - McQuarrie, D. A. (2000). Statistical Mechanics. University Science Books.
# ==============================================================================

@doc """
    sk_to_gr(q_grid::AbstractVector{T}, Sq_grid::AbstractVector{T}, r_grid::AbstractVector{T}, ρ::T) where {T<:AbstractFloat}

Inverts the static structure factor ``S(q)`` to obtain the real-space pair distribution function ``g(r)``:
``g(r) = 1 + \\frac{1}{2\\pi^2 \\rho r} \\int_0^\\infty q [S(q) - 1] \\sin(q r) \\, dq``

# Arguments
- `q_grid::AbstractVector{T}`: Evenly spaced grid of wavevector magnitudes ``q``.
- `Sq_grid::AbstractVector{T}`: Evaluated static structure factor ``S(q)`` at `q_grid`.
- `r_grid::AbstractVector{T}`: Target radial distance grid ``r``.
- `ρ::T`: Number density ``\\rho = N / V``.

# Returns
- `Vector{T}`: Pair distribution function ``g(r)`` evaluated on `r_grid`.
"""
function sk_to_gr(q_grid::AbstractVector{T}, Sq_grid::AbstractVector{T}, r_grid::AbstractVector{T}, ρ::T) where {T<:AbstractFloat}
    if length(q_grid) != length(Sq_grid)
        throw(DimensionMismatch("q_grid and Sq_grid must have the same length."))
    end
    if length(q_grid) < 2
        throw(ArgumentError("q_grid must have at least 2 points."))
    end

    dq = q_grid[2] - q_grid[1]
    gr_result = zeros(T, length(r_grid))
    factor = one(T) / (T(2.0) * (T(π)^2) * ρ)

    # Integrand weights: q * (S(q) - 1)
    f_q = q_grid .* (Sq_grid .- one(T))

    for (j, r) in enumerate(r_grid)
        if r < T(1e-6)
            # Limit r -> 0: (1 / (2*pi^2 * rho)) * \int q^2 (S(q) - 1) dq
            int_val = sum((q_grid .^ 2) .* (Sq_grid .- one(T))) * dq
            gr_result[j] = one(T) + factor * int_val
        else
            int_val = sum(f_q .* sin.(q_grid .* r)) * dq
            gr_result[j] = one(T) + (factor / r) * int_val
        end
    end

    return gr_result
end

@doc """
    gr_to_sk(r_grid::AbstractVector{T}, gr_grid::AbstractVector{T}, q_grid::AbstractVector{T}, ρ::T) where {T<:AbstractFloat}

Calculates the static structure factor ``S(q)`` from the radial distribution function ``g(r)``:
``S(q) = 1 + \\frac{4\\pi \\rho}{q} \\int_0^\\infty r [g(r) - 1] \\sin(q r) \\, dr``

# Arguments
- `r_grid::AbstractVector{T}`: Grid of radial distances ``r``.
- `gr_grid::AbstractVector{T}`: Evaluated pair distribution function ``g(r)``.
- `q_grid::AbstractVector{T}`: Target wavevector grid ``q``.
- `ρ::T`: Number density ``\\rho``.

# Returns
- `Vector{T}`: Static structure factor ``S(q)`` on `q_grid`.
"""
function gr_to_sk(r_grid::AbstractVector{T}, gr_grid::AbstractVector{T}, q_grid::AbstractVector{T}, ρ::T) where {T<:AbstractFloat}
    if length(r_grid) != length(gr_grid)
        throw(DimensionMismatch("r_grid and gr_grid must have the same length."))
    end
    if length(r_grid) < 2
        throw(ArgumentError("r_grid must have at least 2 points."))
    end

    dr = r_grid[2] - r_grid[1]
    Sq_result = zeros(T, length(q_grid))
    factor = T(4.0) * T(π) * ρ

    # Integrand weights: r * (g(r) - 1)
    h_r = r_grid .* (gr_grid .- one(T))

    for (j, q) in enumerate(q_grid)
        if q < T(1e-6)
            # Limit q -> 0: 1 + 4*pi*rho * \int r^2 (g(r) - 1) dr
            int_val = sum((r_grid .^ 2) .* (gr_grid .- one(T))) * dr
            Sq_result[j] = one(T) + factor * int_val
        else
            int_val = sum(h_r .* sin.(q .* r_grid)) * dr
            Sq_result[j] = one(T) + (factor / q) * int_val
        end
    end

    return Sq_result
end

# src/models/sticky_hardsphere.jl
#
# Analytical Structure Factor for Sticky Hard Spheres (SHS)
# Under the Percus-Yevick closure with Wiener-Hopf factorization
#
# References:
# [1] R. J. Baxter, "Percus–Yevick Equation for Hard Spheres with Surface Adhesion",
#     The Journal of Chemical Physics 49(6), 2770–2774 (1968).
#     DOI: 10.1063/1.1670482
# [2] S. V. G. Menon, C. Manohar, and K. S. Rao, "A new interpretation of the sticky hard sphere model",
#     The Journal of Chemical Physics 95(12), 9186–9190 (1991).
#     DOI: 10.1063/1.461199

"""
    baxter_shs_critical_tau(η::Real)

Calculates the critical stickiness parameter τ_c for a given volume fraction η,
corresponding to the fluid-fluid spinodal boundary (where the discriminant D = 0).
The spinodal curve is given by:
    τ_s(η) = (√(η * (1 + η/2) / 3) - η) / (1 - η)
with maximum at the critical point η_c = 3/√2 - 2 ≈ 0.12132, τ_c = (2 - √2)/6 ≈ 0.09763.
"""
function baxter_shs_critical_tau(η::Real)
    T = promote_type(typeof(η), Float64)
    eta = T(η)
    if eta <= zero(T) || eta >= one(T)
        throw(DomainError(η, "Volume fraction η must be in (0, 1)."))
    end
    # Discriminant D = B^2 - 4*A*C = 0 => B = -2*sqrt(A*C)
    # B = -(tau + eta / (1 - eta)) => tau_s = 2*sqrt(A*C) - eta / (1 - eta)
    # = (sqrt(eta * (1 + eta/2) / 3) - eta) / (1 - eta)
    tau_c = (sqrt(eta * (one(T) + eta / T(2.0)) / T(3.0)) - eta) / (one(T) - eta)
    return tau_c
end

"""
    baxter_shs_lambda(η::Real, τ::Real)

Calculates the Baxter stickiness parameter λ by solving the quadratic equation:
    A * λ² + B * λ + C = 0
where:
    A = η / 12
    B = -(τ + η / (1 - η))
    C = (1 + η / 2) / (1 - η)²

Selects the smaller physically acceptable root λ = (-B - √(B² - 4AC)) / (2A).
Throws an error if the system lies within the phase separation spinodal region (discriminant < 0).

# References
- R. J. Baxter, J. Chem. Phys. 49, 2770–2774 (1968).
- S. V. G. Menon, C. Manohar, and K. S. Rao, J. Chem. Phys. 95, 9186–9190 (1991).
"""
function baxter_shs_lambda(η::Real, τ::Real)
    T = promote_type(typeof(η), typeof(τ), Float64)
    eta = T(η)
    tau = T(τ)

    if eta <= zero(T) || eta >= one(T)
        throw(DomainError(η, "Volume fraction η must be in (0, 1)."))
    end
    if tau <= zero(T)
        throw(DomainError(τ, "Stickiness parameter τ must be positive (τ > 0)."))
    end

    A_quad = eta / T(12.0)
    B_quad = -(tau + (eta / (one(T) - eta)))
    C_quad = (one(T) + eta / T(2.0)) / ((one(T) - eta)^2)

    discriminant = B_quad^2 - T(4.0) * A_quad * C_quad

    if discriminant < zero(T)
        tau_c = baxter_shs_critical_tau(eta)
        error("The system lies within the two-phase spinodal separation region: τ = $τ < τ_c ≈ $tau_c for η = $η.")
    end

    # Smaller physical root
    lambda_val = (-B_quad - sqrt(discriminant)) / (T(2.0) * A_quad)
    return lambda_val
end

"""
    _shs_AB_functions(η::T, lambda_val::T, κ::T) where {T<:Real}

Helper function calculating the Baxter functions A(κ) and B(κ) such that
S(κ)⁻¹ = A(κ)² + B(κ)². Includes high-precision Taylor series for |κ| < 0.05.
"""
function _shs_AB_functions(η::T, lambda_val::T, κ::T) where {T<:Real}
    one_minus_eta = one(T) - η
    mu = lambda_val * η * one_minus_eta
    alpha = (one(T) + T(2.0) * η - mu) / (one_minus_eta^2)
    beta = (-T(3.0) * η + mu) / (T(2.0) * (one_minus_eta^2))

    if abs(κ) < T(0.05)
        k2 = κ * κ
        k4 = k2 * k2
        k6 = k4 * k2

        # Taylor expansions for A(κ)
        term_A1 = alpha * (one(T) / T(3.0) - k2 / T(30.0) + k4 / T(840.0) - k6 / T(45360.0))
        term_A2 = beta * (one(T) / T(2.0) - k2 / T(24.0) + k4 / T(720.0) - k6 / T(40320.0))
        term_A3 = -(lambda_val / T(12.0)) * (one(T) - k2 / T(6.0) + k4 / T(120.0) - k6 / T(5040.0))
        A_k = one(T) + T(12.0) * η * (term_A1 + term_A2 + term_A3)

        # Taylor expansions for B(κ)
        term_B1 = alpha * (κ / T(8.0) - (κ * k2) / T(144.0) + (κ * k4) / T(5760.0))
        term_B2 = beta * (κ / T(6.0) - (κ * k2) / T(120.0) + (κ * k4) / T(5040.0))
        term_B3 = -(lambda_val / T(12.0)) * (κ / T(2.0) - (κ * k2) / T(24.0) + (κ * k4) / T(720.0))
        B_k = T(12.0) * η * (term_B1 + term_B2 + term_B3)
    else
        sin_k, cos_k = sincos(κ)

        # Term A(κ)
        term_A1 = alpha * (sin_k - κ * cos_k) / (κ^3)
        term_A2 = beta * (one(T) - cos_k) / (κ^2)
        term_A3 = -(lambda_val / T(12.0)) * (sin_k / κ)
        A_k = one(T) + T(12.0) * η * (term_A1 + term_A2 + term_A3)

        # Term B(κ)
        term_B1 = alpha * (T(0.5) / κ - sin_k / (κ^2) + (one(T) - cos_k) / (κ^3))
        term_B2 = beta * (one(T) / κ - sin_k / (κ^2))
        term_B3 = -(lambda_val / T(12.0)) * ((one(T) - cos_k) / κ)
        B_k = T(12.0) * η * (term_B1 + term_B2 + term_B3)
    end

    return A_k, B_k
end

"""
    S_Baxter_SHS(η::Real, τ::Real, k::Real)

Calculates the static structure factor S(k) for Baxter's Sticky Hard Sphere (SHS) model
under Percus-Yevick closure.

# Arguments
- `η::Real`: Volume fraction of sticky hard spheres (0 < η < 1).
- `τ::Real`: Baxter stickiness parameter (τ > 0). Larger τ corresponds to weaker attraction; τ → ∞ recovers pure hard spheres.
- `k::Real`: Dimensionless wavevector k = q * a (where a is the sphere diameter).

# Returns
- `Float64`: Static structure factor S(k).

# References
- R. J. Baxter, "Percus–Yevick Equation for Hard Spheres with Surface Adhesion", J. Chem. Phys. 49(6), 2770–2774 (1968).
- S. V. G. Menon, C. Manohar, and K. S. Rao, "A new interpretation of the sticky hard sphere model", J. Chem. Phys. 95(12), 9186–9190 (1991).
"""
function S_Baxter_SHS(η::Real, τ::Real, k::Real)
    T = promote_type(typeof(η), typeof(τ), typeof(k), Float64)
    eta = T(η)
    tau = T(τ)
    kappa = T(k)

    lambda_val = baxter_shs_lambda(eta, tau)
    A_k, B_k = _shs_AB_functions(eta, lambda_val, kappa)
    return one(T) / (A_k^2 + B_k^2)
end

"""
    IS_Baxter_SHS(η::Real, τ::Real, k::Real)

Calculates the inverse structure factor 1 / S(k) = A(k)² + B(k)² for Baxter's Sticky Hard Sphere model.
"""
function IS_Baxter_SHS(η::Real, τ::Real, k::Real)
    T = promote_type(typeof(η), typeof(τ), typeof(k), Float64)
    eta = T(η)
    tau = T(τ)
    kappa = T(k)

    lambda_val = baxter_shs_lambda(eta, tau)
    A_k, B_k = _shs_AB_functions(eta, lambda_val, kappa)
    return A_k^2 + B_k^2
end

"""
    S_Menon_SHS(ϕ::Real, u0_kT::Real, Δ::Real, σ::Real, k::Real)

Calculates the static structure factor S(k) for sticky hard spheres using the physical
square-well mapping derived by Menon, Manohar, and Rao (1991).

# Arguments
- `ϕ::Real`: Physical volume fraction of particles.
- `u0_kT::Real`: Dimensionless attractive well depth u₀ / (k_B * T) (typically negative, e.g., -1.5).
- `Δ::Real`: Width of the attractive square well.
- `σ::Real`: Hard-core particle diameter.
- `k::Real`: Physical scattering wavevector (dimension 1 / length).

# Returns
- `Float64`: Static structure factor S(k).

# References
- S. V. G. Menon, C. Manohar, and K. S. Rao, "A new interpretation of the sticky hard sphere model", J. Chem. Phys. 95(12), 9186–9190 (1991).
"""
function S_Menon_SHS(ϕ::Real, u0_kT::Real, Δ::Real, σ::Real, k::Real)
    T = promote_type(typeof(ϕ), typeof(u0_kT), typeof(Δ), typeof(σ), typeof(k), Float64)
    phi = T(ϕ)
    u0 = T(u0_kT)
    delta = T(Δ)
    sigma = T(σ)
    kv = T(k)

    if delta <= zero(T) || sigma <= zero(T)
        throw(DomainError((Δ, σ), "Well width Δ and diameter σ must be positive."))
    end

    # 1. Potential range and perturbation parameter
    a = sigma + delta
    epsilon = delta / a

    # 2. Rescaled density η and stickiness parameter τ (Eq. in Menon et al. 1991)
    eta = phi / ((one(T) - epsilon)^3)
    tau = (one(T) / (T(12.0) * epsilon)) * exp(u0)

    # 3. Dimensionless wavevector κ = k * a
    kappa = kv * a

    return S_Baxter_SHS(eta, tau, kappa)
end

# Vectorized convenience methods
S_Baxter_SHS(η::Real, τ::Real, k_vec::AbstractVector{<:Real}) = [S_Baxter_SHS(η, τ, k) for k in k_vec]
IS_Baxter_SHS(η::Real, τ::Real, k_vec::AbstractVector{<:Real}) = [IS_Baxter_SHS(η, τ, k) for k in k_vec]
S_Menon_SHS(ϕ::Real, u0_kT::Real, Δ::Real, σ::Real, k_vec::AbstractVector{<:Real}) = [S_Menon_SHS(ϕ, u0_kT, Δ, σ, k) for k in k_vec]

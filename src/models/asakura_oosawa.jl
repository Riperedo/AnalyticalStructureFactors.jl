# src/models/asakura_oosawa.jl
#
# Analytical Structure Factors for Colloid-Polymer Mixtures
# Under the Asakura-Oosawa (AO) Model with Density Functional / Percus-Yevick Closure
#
# Reference:
# M. Schmidt, H. Löwen, J. M. Brader, and R. Evans,
# "Density functional theory for a model colloid-polymer mixture: bulk fluid phases",
# J. Phys.: Condens. Matter 14, 9353–9377 (2002).
# DOI: 10.1088/0953-8984/14/40/310

"""
    _ao_I_minus1(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(a, 1) dx x² j₀(kx) / x = [cos(a*k) - cos(k)] / k²
with stable Taylor expansion for k → 0.
"""
function _ao_I_minus1(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4)
        a2 = a * a
        a4 = a2 * a2
        a6 = a4 * a2
        k2 = k * k
        k4 = k2 * k2
        return T(0.5) * (one(T) - a2) + (a4 - one(T)) * k2 / T(24) + (one(T) - a6) * k4 / T(720)
    else
        return (cos(a * k) - cos(k)) / (k * k)
    end
end

"""
    _ao_J_minus1(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(0, a) dx x² j₀(kx) / x = [1 - cos(a*k)] / k²
with stable Taylor expansion for k → 0.
"""
function _ao_J_minus1(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4) || a < T(1e-4)
        a2 = a * a
        a4 = a2 * a2
        a6 = a4 * a2
        k2 = k * k
        k4 = k2 * k2
        return a2 / T(2) - a4 * k2 / T(24) + a6 * k4 / T(720)
    else
        return (one(T) - cos(a * k)) / (k * k)
    end
end

"""
    _ao_I0(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(a, 1) dx x² j₀(kx) = [-sin(a*k) + a*k*cos(a*k) + sin(k) - k*cos(k)] / k³
with stable Taylor expansion for k → 0.
"""
function _ao_I0(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4)
        a3 = a^3
        a5 = a3 * a * a
        a7 = a5 * a * a
        k2 = k * k
        k4 = k2 * k2
        return (one(T) - a3) / T(3) + (a5 - one(T)) * k2 / T(30) + (one(T) - a7) * k4 / T(840)
    else
        return (-sin(a * k) + a * k * cos(a * k) + sin(k) - k * cos(k)) / (k^3)
    end
end

"""
    _ao_J0(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(0, a) dx x² j₀(kx) = [sin(a*k) - a*k*cos(a*k)] / k³
with stable Taylor expansion for k → 0.
"""
function _ao_J0(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4) || a < T(1e-4)
        a3 = a^3
        a5 = a3 * a * a
        a7 = a5 * a * a
        k2 = k * k
        k4 = k2 * k2
        return a3 / T(3) - a5 * k2 / T(30) + a7 * k4 / T(840)
    else
        return (sin(a * k) - a * k * cos(a * k)) / (k^3)
    end
end

"""
    _ao_I1(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(a, 1) dx x³ j₀(kx)
"""
function _ao_I1(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4)
        a4 = a^4
        a6 = a4 * a * a
        a8 = a6 * a * a
        k2 = k * k
        k4 = k2 * k2
        return (one(T) - a4) / T(4) + (a6 - one(T)) * k2 / T(36) + (one(T) - a8) * k4 / T(960)
    else
        return (((a * k)^2 - T(2)) * cos(a * k) + T(2) * k * (sin(k) - a * sin(a * k)) - (k^2 - T(2)) * cos(k)) / (k^4)
    end
end

"""
    _ao_J1(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(0, a) dx x³ j₀(kx)
"""
function _ao_J1(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4) || a < T(1e-4)
        a4 = a^4
        a6 = a4 * a * a
        a8 = a6 * a * a
        k2 = k * k
        k4 = k2 * k2
        return a4 / T(4) - a6 * k2 / T(36) + a8 * k4 / T(960)
    else
        return ((T(2) - (a * k)^2) * cos(a * k) + T(2) * a * k * sin(a * k) - T(2)) / (k^4)
    end
end

"""
    _ao_I3(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(a, 1) dx x⁵ j₀(kx)
"""
function _ao_I3(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4)
        a6 = a^6
        a8 = a6 * a * a
        a10 = a8 * a * a
        k2 = k * k
        k4 = k2 * k2
        return (one(T) - a6) / T(6) + (a8 - one(T)) * k2 / T(48) + (one(T) - a10) * k4 / T(1200)
    else
        ak = a * k
        num = T(4) * k * (a * (T(6) - ak^2) * sin(ak) + (k^2 - T(6)) * sin(k)) +
              (ak^4 - T(12) * ak^2 + T(24)) * cos(ak) -
              (k^4 - T(12) * k^2 + T(24)) * cos(k)
        return num / (k^6)
    end
end

"""
    _ao_J3(a::T, k::T) where {T<:Real}

Computes the integral:
    ∫(0, a) dx x⁵ j₀(kx)
"""
function _ao_J3(a::T, k::T) where {T<:Real}
    if abs(k) < T(1e-4) || a < T(1e-4)
        a6 = a^6
        a8 = a6 * a * a
        a10 = a8 * a * a
        k2 = k * k
        k4 = k2 * k2
        return a6 / T(6) - a8 * k2 / T(48) + a10 * k4 / T(1200)
    else
        ak = a * k
        num = T(4) * a * k * (ak^2 - T(6)) * sin(ak) - (ak^4 - T(12) * ak^2 + T(24)) * cos(ak) + T(24)
        return num / (k^6)
    end
end

# --- Asakura-Oosawa Free Volume & Spinodal Relations ---

"""
    free_volume_fraction_AO(η_c::Real, q::Real)

Calculates the free volume fraction α(η_c, q) for the Asakura-Oosawa model
from free-volume theory (Eq. 69 in Schmidt et al. 2002):
    α(η_c, q) = (1 - η_c) * exp(-[A(q)*γ + B(q)*γ² + C(q)*γ³])
where γ = η_c / (1 - η_c), q = R_p / R_c (polymer-to-colloid size ratio).
"""
function free_volume_fraction_AO(η_c::Real, q::Real)
    T = promote_type(typeof(η_c), typeof(q), Float64)
    η = T(η_c)
    qv = T(q)
    if η < zero(T) || η >= one(T)
        throw(DomainError(η_c, "Colloid packing fraction η_c must be in [0, 1)."))
    end
    γ = η / (one(T) - η)
    A = qv^3 + T(3) * qv^2 + T(3) * qv
    B = T(3) * qv^3 + T(4.5) * qv^2
    C = T(3) * qv^3
    return (one(T) - η) * exp(-(A * γ + B * γ^2 + C * γ^3))
end

"""
    spinodal_AO(η_c::Real, q::Real)

Calculates the polymer packing fraction η_p on the fluid-fluid spinodal line
for a given colloid packing fraction η_c and size ratio q = R_p / R_c (Eq. 73 in Schmidt et al. 2002).
"""
function spinodal_AO(η_c::Real, q::Real)
    T = promote_type(typeof(η_c), typeof(q), Float64)
    η = T(η_c)
    qv = T(q)
    if η <= zero(T) || η >= one(T)
        throw(DomainError(η_c, "Colloid packing fraction η_c must be in (0, 1)."))
    end
    θ1 = one(T) - η
    θ2 = one(T) + T(2) * η
    num = (θ1^4) * θ2 / η
    den = T(12) * (θ1^3) + T(15) * qv * (θ1^2) * θ2 + T(6) * (qv^2) * θ1 * (θ2^2) + (qv^3) * (θ2^3)
    return num / den
end

"""
    polymer_reservoir_to_system_density_AO(η_p_res::Real, η_c::Real, q::Real)

Converts the polymer reservoir packing fraction η_p^r to the system packing fraction η_p:
    η_p = η_p^r * α(η_c, q)
"""
function polymer_reservoir_to_system_density_AO(η_p_res::Real, η_c::Real, q::Real)
    return η_p_res * free_volume_fraction_AO(η_c, q)
end

"""
    polymer_system_to_reservoir_density_AO(η_p::Real, η_c::Real, q::Real)

Converts the system polymer packing fraction η_p to the reservoir packing fraction η_p^r:
    η_p^r = η_p / α(η_c, q)
"""
function polymer_system_to_reservoir_density_AO(η_p::Real, η_c::Real, q::Real)
    α = free_volume_fraction_AO(η_c, q)
    return η_p / α
end

# --- Direct Correlation Functions in Fourier Space ---

"""
    C_AO_HS(η::Real, k::Real)

Fourier transform of the hard-sphere direct correlation function under Percus-Yevick closure.
"""
function C_AO_HS(η::Real, k::Real)
    T = promote_type(typeof(η), typeof(k), Float64)
    η_val = T(η)
    k_val = T(k)
    α0 = -(one(T) + T(2) * η_val)^2 / ((one(T) - η_val)^4)
    α1 = T(6) * η_val * ((one(T) + T(0.5) * η_val)^2) / ((one(T) - η_val)^4)
    α3 = T(0.5) * η_val * α0
    I0 = _ao_I0(zero(T), k_val)
    I1 = _ao_I1(zero(T), k_val)
    I3 = _ao_I3(zero(T), k_val)
    return T(4) * T(π) * (α0 * I0 + α1 * I1 + α3 * I3)
end

"""
    C_AO_star(η_c::Real, q::Real, k::Real)

Fourier transform of the polymer-induced correction C*(k) to the colloid-colloid direct correlation function
(Eq. 80 in Schmidt et al. 2002).
"""
function C_AO_star(η_c::Real, q::Real, k::Real)
    T = promote_type(typeof(η_c), typeof(q), typeof(k), Float64)
    η = T(η_c)
    qv = T(q)
    kv = T(k)
    γ = η / (one(T) - η)
    A = qv^3 + T(3) * qv^2 + T(3) * qv
    B = T(3) * qv^3 + T(4.5) * qv^2
    C = T(3) * qv^3

    τ0 = T(2) * (one(T) + A + T(2) * A * γ + T(4) * B * γ + T(6) * B * γ^2 + T(3) * C * (γ^2) * (T(3) + T(4) * γ))
    τ1 = -T(3) - T(2) * A + T(2) * B / T(3) - T(6) * A * γ - T(8) * B * γ - T(18) * B * γ^2 - T(2) * C * γ * (-one(T) + T(9) * γ * (one(T) + T(2) * γ))
    τ3 = one(T) + T(2) * A * γ + T(6) * B * γ^2 + T(12) * C * γ^3

    pref = -(T(4) * T(π)) / (((one(T) - η)^2) * T(2) * (qv^3))
    I0 = _ao_I0(zero(T), kv)
    I1 = _ao_I1(zero(T), kv)
    I3 = _ao_I3(zero(T), kv)
    return pref * (τ0 * I0 + τ1 * I1 + τ3 * I3)
end

"""
    C_AO_cp(η_c::Real, q::Real, k::Real)

Fourier transform of the colloid-polymer cross direct correlation function C_cp(k)
(Eq. 83-84 in Schmidt et al. 2002).
"""
function C_AO_cp(η_c::Real, q::Real, k::Real)
    T = promote_type(typeof(η_c), typeof(q), typeof(k), Float64)
    η = T(η_c)
    qv = T(q)
    kv = T(k)
    γ = η / (one(T) - η)
    A = qv^3 + T(3) * qv^2 + T(3) * qv
    B = T(3) * qv^3 + T(4.5) * qv^2
    C = T(3) * qv^3

    r1 = T(0.5) * (one(T) - qv)
    r2 = T(0.5) * (one(T) + qv)

    # Core region 1: 0 < r < 0.5*(1 - q)
    C1 = -T(4) * T(π) * (one(T) + γ) * (one(T) + A * γ + T(2) * B * γ^2 + T(3) * C * γ^3) * _ao_J0(r1, kv)

    # Overlap shell region 2: 0.5*(1 - q) < r < 0.5*(1 + q)
    τcp_m1 = -(T(3) * γ / T(32)) * ((one(T) - qv)^2) * ((T(3) + qv + T(3) * γ * (one(T) + qv))^2)
    τcp_0 = T(0.5) * (T(2) + γ * (T(7) + T(3) * γ * (T(5) + T(3) * γ) + T(9) * γ * (qv^2) + ((one(T) + T(3) * γ)^2) * (qv^3) + T(3) * qv * (one(T) + qv)))
    τcp_1 = -(T(3) * γ / T(4)) * (T(5) + T(3) * γ * (T(4) + T(3) * γ) + T(2) * qv + T(6) * γ * qv + (qv + T(3) * γ * qv)^2)
    τcp_3 = T(0.5) * γ * ((one(T) + T(3) * γ)^2)

    diff_Im1 = _ao_I_minus1(r1, kv) - _ao_I_minus1(r2, kv)
    diff_I0 = _ao_I0(r1, kv) - _ao_I0(r2, kv)
    diff_I1 = _ao_I1(r1, kv) - _ao_I1(r2, kv)
    diff_I3 = _ao_I3(r1, kv) - _ao_I3(r2, kv)

    C2 = -T(4) * T(π) * (one(T) + γ) * (τcp_m1 * diff_Im1 + τcp_0 * diff_I0 + τcp_1 * diff_I1 + τcp_3 * diff_I3)

    return C1 + C2
end

"""
    C_AO_cc(η_c::Real, η_p::Real, q::Real, k::Real)

Colloid-colloid direct correlation function in Fourier space:
    C_cc(k) = C_HS(η_c, k) + η_p * C*(η_c, q, k)
"""
function C_AO_cc(η_c::Real, η_p::Real, q::Real, k::Real)
    return C_AO_HS(η_c, k) + η_p * C_AO_star(η_c, q, k)
end

"""
    C_AO_eff_cc(η_c::Real, η_p::Real, q::Real, k::Real)

Effective one-component colloid direct correlation function:
    C_cc^eff(k) = C_cc(k) + ρ_p * [C_cp(k)]²
"""
function C_AO_eff_cc(η_c::Real, η_p::Real, q::Real, k::Real)
    T = promote_type(typeof(η_c), typeof(η_p), typeof(q), typeof(k), Float64)
    ρ_p = T(6) * T(η_p) / (T(π) * (T(q)^3))
    Ccc = C_AO_cc(η_c, η_p, q, k)
    Ccp = C_AO_cp(η_c, q, k)
    return Ccc + ρ_p * (Ccp^2)
end

# --- Partial Structure Factors S_cc, S_cp, S_pp ---

"""
    S_AO_cc(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)

Calculates the Colloid-Colloid partial structure factor S_cc(k) for the Asakura-Oosawa mixture.
Optional `VW = true` applies the Verlet-Weis semi-empirical core correction.
"""
function S_AO_cc(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)
    T = promote_type(typeof(η_c), typeof(η_p), typeof(q), typeof(k), Float64)
    η_in = T(η_c)
    k_in = T(k)
    if VW
        η_in = phi_VW(η_in)
        k_in = k_VW(η_c, k_in)
    end
    ρ_c = T(6) * η_in / T(π)
    ρ_p = T(6) * T(η_p) / (T(π) * (T(q)^3))

    Ccc = C_AO_cc(η_in, η_p, q, k_in)
    Ccp = C_AO_cp(η_in, q, k_in)

    D = one(T) - ρ_c * Ccc - ρ_c * ρ_p * (Ccp^2)
    h_cc = (Ccc + ρ_p * (Ccp^2)) / D
    return one(T) + ρ_c * h_cc
end

"""
    S_AO_cp(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)

Calculates the Colloid-Polymer cross partial structure factor S_cp(k) for the Asakura-Oosawa mixture.
"""
function S_AO_cp(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)
    T = promote_type(typeof(η_c), typeof(η_p), typeof(q), typeof(k), Float64)
    η_in = T(η_c)
    k_in = T(k)
    if VW
        η_in = phi_VW(η_in)
        k_in = k_VW(η_c, k_in)
    end
    ρ_c = T(6) * η_in / T(π)
    ρ_p = T(6) * T(η_p) / (T(π) * (T(q)^3))

    Ccc = C_AO_cc(η_in, η_p, q, k_in)
    Ccp = C_AO_cp(η_in, q, k_in)

    D = one(T) - ρ_c * Ccc - ρ_c * ρ_p * (Ccp^2)
    h_cp = Ccp / D
    return sqrt(ρ_c * ρ_p) * h_cp
end

"""
    S_AO_pp(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)

Calculates the Polymer-Polymer partial structure factor S_pp(k) for the Asakura-Oosawa mixture.
"""
function S_AO_pp(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)
    T = promote_type(typeof(η_c), typeof(η_p), typeof(q), typeof(k), Float64)
    η_in = T(η_c)
    k_in = T(k)
    if VW
        η_in = phi_VW(η_in)
        k_in = k_VW(η_c, k_in)
    end
    ρ_c = T(6) * η_in / T(π)
    ρ_p = T(6) * T(η_p) / (T(π) * (T(q)^3))

    Ccc = C_AO_cc(η_in, η_p, q, k_in)
    Ccp = C_AO_cp(η_in, q, k_in)

    D = one(T) - ρ_c * Ccc - ρ_c * ρ_p * (Ccp^2)
    h_pp = ρ_c * (Ccp^2) / D
    return one(T) + ρ_p * h_pp
end

"""
    S_AO_matrix(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)

Returns the 2x2 Ashcroft-Langreth partial structure factor matrix:
    [ S_cc  S_cp ]
    [ S_cp  S_pp ]
"""
function S_AO_matrix(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)
    s_cc = S_AO_cc(η_c, η_p, q, k; VW = VW)
    s_cp = S_AO_cp(η_c, η_p, q, k; VW = VW)
    s_pp = S_AO_pp(η_c, η_p, q, k; VW = VW)
    return [s_cc s_cp; s_cp s_pp]
end

"""
    S_AO_eff_cc(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)

Calculates the effective one-component colloid structure factor:
    S_cc^eff(k) = 1 / [1 - ρ_c * C_cc^eff(k)]
"""
function S_AO_eff_cc(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)
    T = promote_type(typeof(η_c), typeof(η_p), typeof(q), typeof(k), Float64)
    η_in = T(η_c)
    k_in = T(k)
    if VW
        η_in = phi_VW(η_in)
        k_in = k_VW(η_c, k_in)
    end
    ρ_c = T(6) * η_in / T(π)
    C_eff = C_AO_eff_cc(η_in, η_p, q, k_in)
    return one(T) / (one(T) - ρ_c * C_eff)
end

"""
    S_AO_eff_pp(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)

Calculates the effective polymer structure factor:
    S_pp^eff(k) = 1 / [1 - ρ_p * C_pp^eff(k)]
"""
function S_AO_eff_pp(η_c::Real, η_p::Real, q::Real, k::Real; VW::Bool = false)
    T = promote_type(typeof(η_c), typeof(η_p), typeof(q), typeof(k), Float64)
    η_in = T(η_c)
    k_in = T(k)
    if VW
        η_in = phi_VW(η_in)
        k_in = k_VW(η_c, k_in)
    end
    ρ_c = T(6) * η_in / T(π)
    ρ_p = T(6) * T(η_p) / (T(π) * (T(q)^3))

    Ccc = C_AO_cc(η_in, η_p, q, k_in)
    Ccp = C_AO_cp(η_in, q, k_in)

    C_eff_pp = ρ_c * (Ccp^2) / (one(T) - ρ_c * Ccc)
    return one(T) / (one(T) - ρ_p * C_eff_pp)
end

# Vectorized convenience methods for S(k)
S_AO_cc(η_c::Real, η_p::Real, q::Real, k_vec::AbstractVector{<:Real}; VW::Bool = false) = [S_AO_cc(η_c, η_p, q, k; VW = VW) for k in k_vec]
S_AO_cp(η_c::Real, η_p::Real, q::Real, k_vec::AbstractVector{<:Real}; VW::Bool = false) = [S_AO_cp(η_c, η_p, q, k; VW = VW) for k in k_vec]
S_AO_pp(η_c::Real, η_p::Real, q::Real, k_vec::AbstractVector{<:Real}; VW::Bool = false) = [S_AO_pp(η_c, η_p, q, k; VW = VW) for k in k_vec]
S_AO_eff_cc(η_c::Real, η_p::Real, q::Real, k_vec::AbstractVector{<:Real}; VW::Bool = false) = [S_AO_eff_cc(η_c, η_p, q, k; VW = VW) for k in k_vec]
S_AO_eff_pp(η_c::Real, η_p::Real, q::Real, k_vec::AbstractVector{<:Real}; VW::Bool = false) = [S_AO_eff_pp(η_c, η_p, q, k; VW = VW) for k in k_vec]
S_AO_matrix(η_c::Real, η_p::Real, q::Real, k_vec::AbstractVector{<:Real}; VW::Bool = false) = [S_AO_matrix(η_c, η_p, q, k; VW = VW) for k in k_vec]

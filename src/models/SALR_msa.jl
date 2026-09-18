# src/models/SALR_msa.jl

# ==============================================================================
# EXACT ANALYTICAL TWO-YUKAWA / SALR MEAN SPHERICAL APPROXIMATION (MSA)
# ==============================================================================
# Exact solution of the Ornstein-Zernike equation with the Mean Spherical
# Approximation for fluids interacting via a Hard-Core + Two-Yukawa (SALR) potential:
#   V(r)/(kB*T) = -K1 * exp(-z1*(r-1))/r - K2 * exp(-z2*(r-1))/r   for r > 1 (sigma = 1)
#
# References:
# - Blum, L. (1980). Primitive conditions for the factorizable cluster expansion. J. Stat. Phys., 22(6), 661-672.
# - Cummings, P. T., & Smith, W. R. (1979). Mol. Phys., 38(3), 997-1001.
# - Liu, Y., Chen, W.-R., & Chen, S.-H. (2005). Cluster formation in two-Yukawa fluids. J. Chem. Phys., 122(4), 044507.
# ==============================================================================

using LinearAlgebra: eigvals

# --- Helper functions for direct correlation function C(k) ---

function _tysq_multiTerm(v::T, z::T, bigK::T, k::T) where {T<:AbstractFloat}
    if abs(k) < T(1e-4)
        term1 = (v / z) * (T(0.5) - (one(T) - exp(-z) * (one(T) + z)) / (z^2))
        term2 = (v^2 / (T(4.0) * bigK * (z^2) * exp(z))) * ((T(2.0) - exp(z) * (one(T) - z) - exp(-z) * (one(T) + z)) / (z^2) - one(T))
        term3 = -bigK * (one(T) + z) / (z^2)
        return term1 + term2 + term3
    end
    k2 = k^2; z2 = z^2; ez = exp(z); emz = exp(-z)
    sink = sin(k); cosk = cos(k); zk_sin = (z / k) * sink
    t1 = (v / z) * ((one(T) - cosk) / k2 - (one(T) - emz * (zk_sin + cosk)) / (z2 + k2))
    t2 = (v^2 / (T(4.0) * bigK * z2 * ez)) * ((T(2.0) - ez * (-zk_sin + cosk) - emz * (zk_sin + cosk)) / (z2 + k2) - T(2.0) * (one(T) - cosk) / k2)
    t3 = -(bigK / (z2 + k2)) * (zk_sin + cosk)
    return t1 + t2 + t3
end

function _tysq_Ck(a00::T, b00::T, v1::T, v2::T, k::T, ϕ::T, z1::T, K1::T, z2::T, K2::T) where {T<:AbstractFloat}
    ρ = T(6.0) * ϕ / T(π)
    if abs(k) < T(1e-4)
        return (one(T) - a00) / ρ
    end
    k2 = k^2; k3 = k^3; k4 = k^4; k6 = k^6
    sink = sin(k); cosk = cos(k)
    term_a0 = a00 * (-k * cosk + sink) / k3
    term_b0 = (b00 / k4) * (-k2 * cosk + T(2.0) * k * sink + T(2.0) * (cosk - one(T)))
    term_c0 = (T(0.5) * ϕ * a00 / k6) * (-k4 * cosk + T(4.0) * k3 * sink + T(12.0) * k2 * cosk - T(24.0) * k * sink - T(24.0) * (cosk - one(T)))
    m1 = (abs(K1) > eps(T)) ? _tysq_multiTerm(v1, z1, K1, k) : zero(T)
    m2 = (abs(K2) > eps(T)) ? _tysq_multiTerm(v2, z2, K2, k) : zero(T)
    data = -T(24.0) * ϕ * (term_a0 + term_b0 + term_c0 + m1 + m2)
    return data / ρ
end

# --- Two-Yukawa Baxter Matrix and MSA System Solver ---

function _tysq_calc_coefficients(ϕ::T, K1::T, z1::T, K2::T, z2::T) where {T<:AbstractFloat}
    k1 = K1 * exp(z1)
    k2 = K2 * exp(z2)
    z = [z1, z2]
    k = [k1, k2]

    a0 = (1 + 2*ϕ)/(1 - 2*ϕ + ϕ^2)
    b0 = -3*ϕ/2/(1 - 2*ϕ + ϕ^2)

    v_func(i) = -12*ϕ*(1+2*ϕ)*(1-exp(-z[i])*(1+z[i]))/z[i]/(1-2*ϕ+ϕ^2) +
                8*(1/z[i]^2 - exp(-z[i])*(0.5 + (1+z[i])/z[i]^2)) +
                (8*(1+2*ϕ)*(-1+4*ϕ)*(1/z[i]^2 - exp(-z[i])*(0.5+(1+z[i])/z[i]^2)))/(1-2*ϕ+ϕ^2)
    w_func(i) = 8/z[i]^2 + 8*(1 + 2*ϕ)*(-1 + 4*ϕ)/(1-2*ϕ+ϕ^2)/z[i]^2 - 12*ϕ*(1+2*ϕ)/(1-2*ϕ+ϕ^2)/z[i]
    x_func(i) = (18*ϕ^2*(1-exp(-z[i])*(1+z[i]))/z[i] - (12*ϕ*(-1+4*ϕ)*(1/z[i]^2-exp(-z[i])*(0.5+(1+z[i])/z[i]^2)))) / (1-2*ϕ+ϕ^2)
    y_func(i) = -12*ϕ*(-1+4*ϕ)/(1-2*ϕ+ϕ^2)/z[i]^2 + 18*ϕ^2/(1-2*ϕ+ϕ^2)/z[i]

    Ccd1_11 = -6*ϕ - 6*exp(-2*z[1])*ϕ + 12*exp(-z[1])*ϕ + 6*ϕ*v_func(1) + 12*ϕ*x_func(1) - 12*ϕ*v_func(1)/z[1]^2 + 12*exp(-z[1])*ϕ*v_func(1)/z[1]^2 + 12*exp(-z[1])*ϕ*v_func(1)/z[1] - 12*ϕ*x_func(1)/z[1] + 12*exp(-z[1])*ϕ*x_func(1)/z[1]
    Ccd1_21 = 12*exp(-z[2])*ϕ + 6*ϕ*v_func(2) + 12*ϕ*x_func(2) - 12*ϕ*v_func(2)/z[1]^2 + 12*exp(-z[1])*ϕ*v_func(2)/z[1]^2 + 12*exp(-z[1])*ϕ*v_func(2)/z[1] - 12*ϕ*x_func(2)/z[1] + 12*exp(-z[1])*ϕ*x_func(2)/z[1] - 12*ϕ*z[1]/(z[1]+z[2]) - 12*exp(-z[2])*exp(-z[1])*ϕ*z[2]/(z[1]+z[2])
    Ccd2_12 = 12*exp(-z[1])*ϕ + 6*ϕ*v_func(1) + 12*ϕ*x_func(1) - 12*ϕ*v_func(1)/z[2]^2 + 12*exp(-z[2])*ϕ*v_func(1)/z[2]^2 + 12*exp(-z[2])*ϕ*v_func(1)/z[2] - 12*ϕ*x_func(1)/z[2] + 12*exp(-z[2])*ϕ*x_func(1)/z[2] - 12*ϕ*z[2]/(z[1]+z[2]) - 12*exp(-z[1]-z[2])*ϕ*z[1]/(z[1]+z[2])
    Ccd2_22 = -6*ϕ + 12*exp(-z[2])*ϕ + 6*ϕ*v_func(2) + 12*ϕ*x_func(2) - 12*ϕ*v_func(2)/z[2]^2 + 12*exp(-z[2])*ϕ*v_func(2)/z[2]^2 + 12*exp(-z[2])*ϕ*v_func(2)/z[2] - 12*ϕ*x_func(2)/z[2] + 12*exp(-z[2])*ϕ*x_func(2)/z[2] - 6*exp(-z[2])*exp(-z[2])*ϕ

    Cdd1_11 = -6*ϕ + 6*ϕ*w_func(1) + 12*ϕ*y_func(1) - 12*ϕ*w_func(1)/z[1]^2 + 12*exp(-z[1])*ϕ*w_func(1)/z[1]^2 + 12*exp(-z[1])*ϕ*w_func(1)/z[1] - 12*ϕ*y_func(1)/z[1] + 12*exp(-z[1])*ϕ*y_func(1)/z[1]
    Cdd1_12 = 6*ϕ*w_func(2) + 12*ϕ*y_func(2) - 12*ϕ*w_func(2)/z[1]^2 + 12*exp(-z[1])*ϕ*w_func(2)/z[1]^2 + 12*exp(-z[1])*ϕ*w_func(2)/z[1] - 12*ϕ*y_func(2)/z[1] + 12*exp(-z[1])*ϕ*y_func(2)/z[1] - 12*ϕ*z[1]/(z[1]+z[2])
    Cdd2_12 = 6*ϕ*w_func(1) + 12*ϕ*y_func(1) - 12*ϕ*w_func(1)/z[2]^2 + 12*exp(-z[2])*ϕ*w_func(1)/z[2]^2 + 12*exp(-z[2])*ϕ*w_func(1)/z[2] - 12*ϕ*y_func(1)/z[2] + 12*exp(-z[2])*ϕ*y_func(1)/z[2] - 12*ϕ*z[2]/(z[1]+z[2])
    Cdd2_22 = -6*ϕ + 6*ϕ*w_func(2) + 12*ϕ*y_func(2) - 12*ϕ*w_func(2)/z[2]^2 + 12*exp(-z[2])*ϕ*w_func(2)/z[2]^2 + 12*exp(-z[2])*ϕ*w_func(2)/z[2] - 12*ϕ*y_func(2)/z[2] + 12*exp(-z[2])*ϕ*y_func(2)/z[2]

    Cd1_1 = 6*a0*ϕ + 12*b0*ϕ - 12*a0*ϕ/z[1]^2 + 12*a0*exp(-z[1])*ϕ/z[1]^2 - 12*b0*ϕ/z[1] + 12*a0*exp(-z[1])*ϕ/z[1] + 12*b0*exp(-z[1])*ϕ/z[1] + z[1]
    Cd2_2 = 6*a0*ϕ + 12*b0*ϕ - 12*a0*ϕ/z[2]^2 + 12*a0*exp(-z[2])*ϕ/z[2]^2 - 12*b0*ϕ/z[2] + 12*a0*exp(-z[2])*ϕ/z[2] + 12*b0*exp(-z[2])*ϕ/z[2] + z[2]

    aFNumd01 = -Ccd2_22*k[1]*v_func(1) + Ccd2_12*k[1]*v_func(2)
    aFNumd10 = Ccd1_21*k[2]*v_func(1) - Ccd1_11*k[2]*v_func(2)
    aFNumd11 = a0*Ccd1_21*Ccd2_12 - a0*Ccd1_11*Ccd2_22 + Ccd2_22*Cd1_1*v_func(1) - Ccd1_21*Cd2_2*v_func(1) - Ccd2_12*Cd1_1*v_func(2) + Ccd1_11*Cd2_2*v_func(2)
    aFNumd12 = Ccd2_22*Cdd1_12*v_func(1) - Ccd1_21*Cdd2_22*v_func(1) - Ccd2_12*Cdd1_12*v_func(2) + Ccd1_11*Cdd2_22*v_func(2) + Ccd1_21*Ccd2_12*w_func(2) - Ccd1_11*Ccd2_22*w_func(2)
    aFNumd21 = Ccd2_22*Cdd1_11*v_func(1) - Ccd1_21*Cdd2_12*v_func(1) - Ccd2_12*Cdd1_11*v_func(2) + Ccd1_11*Cdd2_12*v_func(2) + Ccd1_21*Ccd2_12*w_func(1) - Ccd1_11*Ccd2_22*w_func(1)

    bFNumd10 = Ccd1_21*k[2]*x_func(1) - Ccd1_11*k[2]*x_func(2)
    bFNumd11 = b0*Ccd1_21*Ccd2_12 - b0*Ccd1_11*Ccd2_22 + Ccd2_22*Cd1_1*x_func(1) - Ccd1_21*Cd2_2*x_func(1) - Ccd2_12*Cd1_1*x_func(2) + Ccd1_11*Cd2_2*x_func(2)
    bFNumd12 = Ccd2_22*Cdd1_12*x_func(1) - Ccd1_21*Cdd2_22*x_func(1) - Ccd2_12*Cdd1_12*x_func(2) + Ccd1_11*Cdd2_22*x_func(2) + Ccd1_21*Ccd2_12*y_func(2) - Ccd1_11*Ccd2_22*y_func(2)
    bFNumd21 = Ccd2_22*Cdd1_11*x_func(1) - Ccd1_21*Cdd2_12*x_func(1) - Ccd2_12*Cdd1_11*x_func(2) + Ccd1_11*Cdd2_12*x_func(2) + Ccd1_21*Ccd2_12*y_func(1) - Ccd1_11*Ccd2_22*y_func(1)

    bNumd01 = -Ccd2_22*k[1]*x_func(1) + Ccd2_12*k[1]*x_func(2)
    abFDend11 = Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22
    c12FD = -abFDend11

    c1F01 = Ccd2_22*k[1]; c1F10 = -Ccd1_21*k[2]; c1F11 = -Ccd2_22*Cd1_1 + Ccd1_21*Cd2_2
    c1F12 = -Ccd2_22*Cdd1_12 + Ccd1_21*Cdd2_22; c1F21 = -Ccd2_22*Cdd1_11 + Ccd1_21*Cdd2_12
    c2F01 = -Ccd2_12*k[1]; c2F10 = Ccd1_11*k[2]; c2F11 = Ccd2_12*Cd1_1 - Ccd1_11*Cd2_2
    c2F12 = Ccd2_12*Cdd1_12 - Ccd1_11*Cdd2_22; c2F21 = Ccd2_12*Cdd1_11 - Ccd1_11*Cdd2_12

    sigma_d01(s) = (aFNumd01/(abFDend11*s^3) + bNumd01/(abFDend11*s^2) - aFNumd01/(2*abFDend11*s) - bNumd01/(abFDend11*s) - (exp(-z[2])*Ccd2_12*k[1])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (exp(-z[1])*Ccd2_22*k[1])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (Ccd2_22*k[1])/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) + (Ccd2_12*k[1])/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])))
    sigma_d10(s) = (aFNumd10/(abFDend11*s^3) + bFNumd10/(abFDend11*s^2) - aFNumd10/(2*abFDend11*s) - bFNumd10/(abFDend11*s) + (exp(-z[2])*Ccd1_11*k[2])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[1])*Ccd1_21*k[2])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (Ccd1_21*k[2])/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd1_11*k[2])/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])))
    sigma_d11(s) = (aFNumd11/(abFDend11*s^3) + bFNumd11/(abFDend11*s^2) - aFNumd11/(2*abFDend11*s) - bFNumd11/(abFDend11*s) + (exp(-z[2])*Ccd2_12*Cd1_1)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[1])*Ccd2_22*Cd1_1)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[2])*Ccd1_11*Cd2_2)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (exp(-z[1])*Ccd1_21*Cd2_2)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (Ccd2_22*Cd1_1)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd1_21*Cd2_2)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd2_12*Cd1_1)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + (Ccd1_11*Cd2_2)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])))
    sigma_d12(s) = (aFNumd12/(abFDend11*s^3) + bFNumd12/(abFDend11*s^2) - aFNumd12/(2*abFDend11*s) - bFNumd12/(abFDend11*s) + (exp(-z[2])*Ccd2_12*Cdd1_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[1])*Ccd2_22*Cdd1_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[2])*Ccd1_11*Cdd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (exp(-z[1])*Ccd1_21*Cdd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (Ccd2_22*Cdd1_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd1_21*Cdd2_22)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd2_12*Cdd1_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + (Ccd1_11*Cdd2_22)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + 1/(s + z[2]))
    sigma_d21(s) = (aFNumd21/(abFDend11*s^3) + bFNumd21/(abFDend11*s^2) - aFNumd21/(2*abFDend11*s) - bFNumd21/(abFDend11*s) + (exp(-z[2])*Ccd2_12*Cdd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[1])*Ccd2_22*Cdd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) - (exp(-z[2])*Ccd1_11*Cdd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + (exp(-z[1])*Ccd1_21*Cdd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) + 1/(s + z[1]) + (Ccd2_22*Cdd1_11)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd1_21*Cdd2_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (Ccd2_12*Cdd1_11)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + (Ccd1_11*Cdd2_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])))

    tau_d01(s) = (aFNumd01*(1/s^3 + 1/s^2)/abFDend11 + bNumd01/(abFDend11*s^2) - (exp(-z[2])*k[1]*z[2]*Ccd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + (exp(-z[1])*k[1]*z[1]*Ccd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])))
    tau_d10(s) = (aFNumd10*(1/s^3 + 1/s^2)/abFDend11 + bFNumd10/(abFDend11*s^2) + (exp(-z[2])*k[2]*z[2]*Ccd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) - (exp(-z[1])*k[2]*z[1]*Ccd1_21)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])))
    tau_d11(s) = (aFNumd11*(1/s^3 + 1/s^2)/abFDend11 + bFNumd11/(abFDend11*s^2) + (exp(-z[2])*Cd1_1*z[2]*Ccd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) - (exp(-z[1])*Cd1_1*z[1]*Ccd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (exp(-z[2])*Cd2_2*z[2]*Ccd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + (exp(-z[1])*Cd2_2*z[1]*Ccd1_21)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])))
    tau_d12(s) = (aFNumd12*(1/s^3 + 1/s^2)/abFDend11 + bFNumd12/(abFDend11*s^2) + (exp(-z[2])*Cdd1_12*z[2]*Ccd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) - (exp(-z[1])*Cdd1_12*z[1]*Ccd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) - (exp(-z[2])*Cdd2_22*z[2]*Ccd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) + (exp(-z[1])*Cdd2_22*z[1]*Ccd1_21)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])))
    tau_d21(s) = (aFNumd21*(1/s^3 + 1/s^2)/abFDend11 + bFNumd21/(abFDend11*s^2) - (exp(-z[1])*Ccd2_22*Cdd1_11*z[1])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) + (exp(-z[1])*Ccd1_21*Cdd2_12*z[1])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[1])) + (exp(-z[2])*Ccd2_12*Cdd1_11*z[2])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])) - (exp(-z[2])*Ccd1_11*Cdd2_12*z[2])/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z[2])))

    # Exact MATLAB E1d coefficients:
    E1d02 = 12*c1F01*ϕ*sigma_d01(z[1]) - 12*c1F01*exp(-z[1])*ϕ*tau_d01(z[1])
    E1d11 = 12*c1F10*ϕ*sigma_d01(z[1]) + 12*c1F01*ϕ*sigma_d10(z[1]) - 12*c1F10*exp(-z[1])*ϕ*tau_d01(z[1]) - 12*c1F01*exp(-z[1])*ϕ*tau_d10(z[1])
    E1d12 = (-c1F01) + 12*c1F11*ϕ*sigma_d01(z[1]) + 12*c1F01*ϕ*sigma_d11(z[1]) - 12*c1F11*exp(-z[1])*ϕ*tau_d01(z[1]) - 12*c1F01*exp(-z[1])*ϕ*tau_d11(z[1])
    E1d13 = 12*c1F12*ϕ*sigma_d01(z[1]) + 12*c1F01*ϕ*sigma_d12(z[1]) - 12*c1F12*exp(-z[1])*ϕ*tau_d01(z[1]) - 12*c1F01*exp(-z[1])*ϕ*tau_d12(z[1])
    E1d20 = 12*c1F10*ϕ*sigma_d10(z[1]) - 12*c1F10*exp(-z[1])*ϕ*tau_d10(z[1])
    E1d21 = (-c1F10) + 12*c1F11*ϕ*sigma_d10(z[1]) + 12*c1F10*ϕ*sigma_d11(z[1]) - 12*c1F11*exp(-z[1])*ϕ*tau_d10(z[1]) - 12*c1F10*exp(-z[1])*ϕ*tau_d11(z[1])
    E1d22 = (-c1F11) + 12*c12FD*ϕ*sigma_d01(z[1]) + 12*c1F21*ϕ*sigma_d01(z[1]) + 12*c1F12*ϕ*sigma_d10(z[1]) + 12*c1F11*ϕ*sigma_d11(z[1]) + 12*c1F10*ϕ*sigma_d12(z[1]) + 12*c1F01*ϕ*sigma_d21(z[1]) - 12*c1F21*exp(-z[1])*ϕ*tau_d01(z[1]) - 12*c1F12*exp(-z[1])*ϕ*tau_d10(z[1]) - 12*c1F11*exp(-z[1])*ϕ*tau_d11(z[1]) - 12*c1F10*exp(-z[1])*ϕ*tau_d12(z[1]) - 12*c1F01*exp(-z[1])*ϕ*tau_d21(z[1])
    E1d23 = (-c1F12) + 12*c1F12*ϕ*sigma_d11(z[1]) + 12*c1F11*ϕ*sigma_d12(z[1]) - 12*c1F12*exp(-z[1])*ϕ*tau_d11(z[1]) - 12*c1F11*exp(-z[1])*ϕ*tau_d12(z[1])
    E1d24 = 12*c1F12*ϕ*sigma_d12(z[1]) - 12*c1F12*exp(-z[1])*ϕ*tau_d12(z[1])
    E1d31 = 12*c12FD*ϕ*sigma_d10(z[1]) + 12*c1F21*ϕ*sigma_d10(z[1]) + 12*c1F10*ϕ*sigma_d21(z[1]) - 12*c1F21*exp(-z[1])*ϕ*tau_d10(z[1]) - 12*c1F10*exp(-z[1])*ϕ*tau_d21(z[1])
    E1d32 = (-c12FD) - c1F21 + 12*c12FD*ϕ*sigma_d11(z[1]) + 12*c1F21*ϕ*sigma_d11(z[1]) + 12*c1F11*ϕ*sigma_d21(z[1]) - 12*c1F21*exp(-z[1])*ϕ*tau_d11(z[1]) - 12*c1F11*exp(-z[1])*ϕ*tau_d21(z[1])
    E1d33 = 12*c12FD*ϕ*sigma_d12(z[1]) + 12*c1F21*ϕ*sigma_d12(z[1]) + 12*c1F12*ϕ*sigma_d21(z[1]) - 12*c1F21*exp(-z[1])*ϕ*tau_d12(z[1]) - 12*c1F12*exp(-z[1])*ϕ*tau_d21(z[1])
    E1d42 = 12*c12FD*ϕ*sigma_d21(z[1]) + 12*c1F21*ϕ*sigma_d21(z[1]) - 12*c1F21*exp(-z[1])*ϕ*tau_d21(z[1])

    # Exact MATLAB E2d coefficients:
    E2d02 = 12*c2F01*exp(z[2])*ϕ*sigma_d01(z[2]) - 12*c2F01*ϕ*tau_d01(z[2])
    E2d11 = 12*c2F10*exp(z[2])*ϕ*sigma_d01(z[2]) + 12*c2F01*exp(z[2])*ϕ*sigma_d10(z[2]) - 12*c2F10*ϕ*tau_d01(z[2]) - 12*c2F01*ϕ*tau_d10(z[2])
    E2d12 = (-c2F01)*exp(z[2]) + 12*c2F11*exp(z[2])*ϕ*sigma_d01(z[2]) + 12*c2F01*exp(z[2])*ϕ*sigma_d11(z[2]) - 12*c2F11*ϕ*tau_d01(z[2]) - 12*c2F01*ϕ*tau_d11(z[2])
    E2d13 = 12*c12FD*exp(z[2])*ϕ*sigma_d01(z[2]) + 12*c2F12*exp(z[2])*ϕ*sigma_d01(z[2]) + 12*c2F01*exp(z[2])*ϕ*sigma_d12(z[2]) - 12*c2F12*ϕ*tau_d01(z[2]) - 12*c2F01*ϕ*tau_d12(z[2])
    E2d20 = 12*c2F10*exp(z[2])*ϕ*sigma_d10(z[2]) - 12*c2F10*ϕ*tau_d10(z[2])
    E2d21 = (-c2F10)*exp(z[2]) + 12*c2F11*exp(z[2])*ϕ*sigma_d10(z[2]) + 12*c2F10*exp(z[2])*ϕ*sigma_d11(z[2]) - 12*c2F11*ϕ*tau_d10(z[2]) - 12*c2F10*ϕ*tau_d11(z[2])
    E2d22 = (-c2F11)*exp(z[2]) + 12*c2F21*exp(z[2])*ϕ*sigma_d01(z[2]) + 12*c12FD*exp(z[2])*ϕ*sigma_d10(z[2]) + 12*c2F12*exp(z[2])*ϕ*sigma_d10(z[2]) + 12*c2F11*exp(z[2])*ϕ*sigma_d11(z[2]) + 12*c2F10*exp(z[2])*ϕ*sigma_d12(z[2]) + 12*c2F01*exp(z[2])*ϕ*sigma_d21(z[2]) - 12*c2F21*ϕ*tau_d01(z[2]) - 12*c2F12*ϕ*tau_d10(z[2]) - 12*c2F11*ϕ*tau_d11(z[2]) - 12*c2F10*ϕ*tau_d12(z[2]) - 12*c2F01*ϕ*tau_d21(z[2])
    E2d23 = (-c12FD)*exp(z[2]) - c2F12*exp(z[2]) + 12*c12FD*exp(z[2])*ϕ*sigma_d11(z[2]) + 12*c2F12*exp(z[2])*ϕ*sigma_d11(z[2]) + 12*c2F11*exp(z[2])*ϕ*sigma_d12(z[2]) - 12*c2F12*ϕ*tau_d11(z[2]) - 12*c2F11*ϕ*tau_d12(z[2])
    E2d24 = 12*c12FD*exp(z[2])*ϕ*sigma_d12(z[2]) + 12*c2F12*exp(z[2])*ϕ*sigma_d12(z[2]) - 12*c2F12*ϕ*tau_d12(z[2])
    E2d31 = 12*c2F21*exp(z[2])*ϕ*sigma_d10(z[2]) + 12*c2F10*exp(z[2])*ϕ*sigma_d21(z[2]) - 12*c2F21*ϕ*tau_d10(z[2]) - 12*c2F10*ϕ*tau_d21(z[2])
    E2d32 = (-c2F21)*exp(z[2]) + 12*c2F21*exp(z[2])*ϕ*sigma_d11(z[2]) + 12*c2F11*exp(z[2])*ϕ*sigma_d21(z[2]) - 12*c2F21*ϕ*tau_d11(z[2]) - 12*c2F11*ϕ*tau_d21(z[2])
    E2d33 = 12*c2F21*exp(z[2])*ϕ*sigma_d12(z[2]) + 12*c12FD*exp(z[2])*ϕ*sigma_d21(z[2]) + 12*c2F12*exp(z[2])*ϕ*sigma_d21(z[2]) - 12*c2F21*ϕ*tau_d12(z[2]) - 12*c2F12*ϕ*tau_d21(z[2])
    E2d42 = 12*c2F21*exp(z[2])*ϕ*sigma_d21(z[2]) - 12*c2F21*ϕ*tau_d21(z[2])

    Ecoeff = T[
        E1d20, E1d11, E1d21, E1d31, E1d02, E1d12, E1d42, E1d13, E1d23, E1d33, E1d24,
        E2d20, E2d11, E2d21, E2d31, E2d02, E2d12, E2d42, E2d13, E2d23, E2d33, E2d24,
        E1d22, E1d32, E2d22, E2d32
    ]

    function ABC12(d1, d2)
        A = ((aFNumd01*d2 + aFNumd21*d1^2*d2 + aFNumd10*d1 + aFNumd11*d1*d2 + aFNumd12*d1*d2^2)/(abFDend11*d1*d2))
        B = ((bNumd01*d2 + bFNumd21*d1^2*d2 + bFNumd10*d1 + bFNumd11*d1*d2 + bFNumd12*d1*d2^2)/(abFDend11*d1*d2))
        C1 = ((-((Ccd2_22*Cd1_1*d1*d2 - Ccd1_21*Cd2_2*d1*d2 + Ccd2_22*Cdd1_11*d1^2*d2 - Ccd1_21*Cdd2_12*d1^2*d2 + Ccd2_22*Cdd1_12*d1*d2^2 - Ccd1_21*Cdd2_22*d1*d2^2 - Ccd2_22*d2*k[1] + Ccd1_21*d1*k[2]))) / (d1*((-Ccd1_21)*Ccd2_12*d2 + Ccd1_11*Ccd2_22*d2)))
        C2 = ((-((Ccd2_12*d2*(((-Cd1_1)*d1 - Cdd1_11*d1^2 - Cdd1_12*d1*d2 + k[1])) - Ccd1_11*d1*(((-Cd2_2)*d2 - Cdd2_12*d1*d2 - Cdd2_22*d2^2 + k[2]))))) / ((-Ccd1_21)*Ccd2_12*d1*d2 + Ccd1_11*Ccd2_22*d1*d2))
        return A, B, C1, C2
    end

    function gHat(a, b, c1, d1, c2, d2, s)
        sigma = a/s^3 + b/s^2 - (a/2 + b + c1*exp(-z[1]) + c2*exp(-z[2]))/s + (c1+d1)/(z[1]+s) + (c2+d2)/(z[2]+s)
        tau = a*(1/s^3 + 1/s^2) + b/s^2 - (z[1]*c1*exp(-z[1])/(z[1]+s) + z[2]*c2*exp(-z[2])/(z[2]+s))/s
        q = sigma - tau*exp(-s)
        return ((a*(s+1) + b*s)/s^2 - z[1]*c1*exp(-z[1])/(z[1]+s) - z[2]*c2*exp(-z[2])/(z[2]+s))*exp(-s) / (1 - 12*ϕ*q)
    end

    function CxCoef(d1, d2)
        a, b, c1, c2 = ABC12(d1, d2)
        ghat1 = gHat(a, b, c1, d1, c2, d2, z[1])
        ghat2 = gHat(a, b, c1, d1, c2, d2, z[2])
        v1_c = 24 * ϕ * K1 * exp(z[1]) * ghat1
        v2_c = 24 * ϕ * K2 * exp(z[2]) * ghat2
        a00 = a^2
        b00 = -12 * ϕ * ((a + b)^2 / 2 + a * (c1 * exp(-z[1]) + c2 * exp(-z[2])))
        return (a=a, b=b, c1=c1, c2=c2, v1=v1_c, v2=v2_c, a00=a00, b00=b00)
    end

    return (Ecoeff=Ecoeff, CxCoef=CxCoef, z=z, k=k)
end

function _tysq_solve_physical_coefficients(ϕ::T, K1::T, z1::T, K2::T, z2::T; allow_swap::Bool=true) where {T<:AbstractFloat}
    # Single Yukawa / Hard-Sphere limits
    if abs(K1) < eps(T) && abs(K2) < eps(T)
        a = (one(T) + T(2.0)*ϕ) / (one(T) - ϕ)^2
        b = -T(1.5)*ϕ / (one(T) - ϕ)^2
        return (a00=a^2, b00=-T(12.0)*ϕ*(a+b)^2/T(2.0), v1=zero(T), v2=zero(T), z1=z1, K1=zero(T), z2=z2, K2=zero(T))
    elseif abs(K1) < T(1e-6)
        c2 = _yukawa_msa_coefficients(ϕ, K2, z2)
        v2 = -T(12.0) * ϕ * c2.c
        a00 = c2.a^2
        b00 = -T(12.0) * ϕ * ((c2.a + c2.b)^2 / T(2.0) + c2.a * c2.c * exp(-z2))
        return (a00=a00, b00=b00, v1=zero(T), v2=v2, z1=z1, K1=zero(T), z2=z2, K2=K2)
    elseif abs(K2) < T(1e-6)
        c1 = _yukawa_msa_coefficients(ϕ, K1, z1)
        v1 = -T(12.0) * ϕ * c1.c
        a00 = c1.a^2
        b00 = -T(12.0) * ϕ * ((c1.a + c1.b)^2 / T(2.0) + c1.a * c1.c * exp(-z1))
        return (a00=a00, b00=b00, v1=v1, v2=zero(T), z1=z1, K1=K1, z2=z2, K2=zero(T))
    elseif abs(z1 - z2) < T(1e-3)
        K_eff = K1 + K2
        c = _yukawa_msa_coefficients(ϕ, K_eff, z1)
        v = -T(12.0) * ϕ * c.c
        a00 = c.a^2
        b00 = -T(12.0) * ϕ * ((c.a + c.b)^2 / T(2.0) + c.a * c.c * exp(-z1))
        return (a00=a00, b00=b00, v1=v, v2=zero(T), z1=z1, K1=K_eff, z2=z2, K2=zero(T))
    end

    calc = _tysq_calc_coefficients(ϕ, K1, z1, K2, z2)
    E = copy(calc.Ecoeff)
    # Zero out unphysical terms exactly as in Blum / Liu / MATLAB
    E[1] = zero(T); E[3] = zero(T); E[9] = zero(T); E[11] = zero(T)
    E[16] = zero(T); E[17] = zero(T); E[26] = zero(T); E[18] = zero(T)

    gE1d42 = E[7]
    gE1d31, gE1d32, gE1d33 = E[4], E[24], E[10]
    gE1d22 = E[23]
    gE1d11, gE1d12, gE1d13 = E[2], E[6], E[8]
    gE1d02 = E[5]

    gE2d20, gE2d11, gE2d21, gE2d31 = E[12], E[13], E[14], E[15]
    gE2d13, gE2d23, gE2d33, gE2d24, gE2d22 = E[19], E[20], E[21], E[22], E[25]

    function eval_res(d2::T, sgn::T)
        abs(d2) < eps(T) && return (valid=false, res=T(NaN), d1=T(NaN))
        y22 = gE2d31 * d2 + gE2d33 * d2^3
        y21 = gE2d20 + gE2d21 * d2 + gE2d22 * d2^2 + gE2d23 * d2^3 + gE2d24 * d2^4
        y20 = gE2d11 * d2 + gE2d13 * d2^3
        disc = y21^2 - T(4.0) * y22 * y20
        disc < zero(T) && return (valid=false, res=T(NaN), d1=T(NaN))
        abs(y22) < eps(T) && return (valid=false, res=T(NaN), d1=T(NaN))
        d1 = (-y21 + sgn * sqrt(disc)) / (T(2.0) * y22)
        abs(d1) < eps(T) && return (valid=false, res=T(NaN), d1=T(NaN))

        y14 = gE1d42 * d2
        y13 = gE1d31 + gE1d32 * d2 + gE1d33 * d2^2
        y12 = gE1d22 * d2
        y11 = gE1d11 + gE1d12 * d2 + gE1d13 * d2^2
        y10 = gE1d02 * d2

        res = y14 * d1^2 + y13 * d1 + y12 + y11 / d1 + y10 / d1^2
        return (valid=true, res=res, d1=d1)
    end

    q_grid = range(T(0.01), T(60.0), length=200)
    r_core = range(T(0.1), T(0.9), length=15)
    ρ = T(6.0) * ϕ / T(π)
    dq = step(q_grid)

    candidates = []
    max_d2 = max(T(50.0), T(10.0) * (abs(K1)*exp(z1) + abs(K2)*exp(z2)))
    min_d2 = T(1e-5)
    log_pts = exp.(range(log(min_d2), log(max_d2), length=400))
    grids = (-reverse(log_pts), log_pts)

    for sgn in (one(T), -one(T))
        for g in grids
            for i in 1:(length(g)-1)
                d2_a, d2_b = g[i], g[i+1]
                va = eval_res(d2_a, sgn)
                vb = eval_res(d2_b, sgn)
                if va.valid && vb.valid && (va.res * vb.res <= zero(T))
                    a, b = d2_a, d2_b
                    fa = va.res
                    for _ in 1:45
                        mid = T(0.5) * (a + b)
                        vm = eval_res(mid, sgn)
                        if !vm.valid || abs(vm.res) < T(1e-14)
                            a = mid; b = mid; break
                        end
                        if fa * vm.res <= zero(T)
                            b = mid
                        else
                            a = mid
                            fa = vm.res
                        end
                    end
                    root_d2 = T(0.5) * (a + b)
                    sol = eval_res(root_d2, sgn)
                    if sol.valid && abs(sol.res) < T(1e-3)
                        coe = calc.CxCoef(sol.d1, root_d2)
                        if coe.a > zero(T) && (coe.v1 * K1 >= -T(1e-5)) && (coe.v2 * K2 >= -T(1e-5))
                            hk_vals = T[]
                            for q in q_grid
                                Ck_val = _tysq_Ck(coe.a00, coe.b00, coe.v1, coe.v2, q, ϕ, z1, K1, z2, K2)
                                sq = one(T) / (one(T) - ρ * Ck_val)
                                push!(hk_vals, (sq - one(T)) / ρ)
                            end
                            core_penalty = zero(T)
                            for r_val in r_core
                                int_val = zero(T)
                                for (iq, q) in enumerate(q_grid)
                                    int_val += hk_vals[iq] * q * sin(q * r_val) * dq
                                end
                                gr = one(T) + int_val / (T(2.0) * T(π)^2 * r_val)
                                core_penalty += abs(gr)
                            end
                            core_penalty /= length(r_core)
                            push!(candidates, (coe=coe, penalty=core_penalty, res1=sol.res))
                        end
                    end
                end
            end
        end
    end

    if isempty(candidates)
        if allow_swap
            try
                swapped = _tysq_solve_physical_coefficients(ϕ, K2, z2, K1, z1; allow_swap=false)
                return (a00=swapped.a00, b00=swapped.b00, v1=swapped.v2, v2=swapped.v1, z1=z1, K1=K1, z2=z2, K2=K2)
            catch
            end
        end
        throw(DomainError((ϕ, K1, z1, K2, z2), "No physical MSA roots found."))
    end

    sort!(candidates, by = c -> c.penalty)
    best = candidates[1]
    return (a00=best.coe.a00, b00=best.coe.b00, v1=best.coe.v1, v2=best.coe.v2, z1=z1, K1=K1, z2=z2, K2=K2)
end

# ------------------------------------------------------------------------------
# PUBLIC API: S_SALR_MSA, IS_SALR_MSA, C_SALR_MSA
# ------------------------------------------------------------------------------

"""
    S_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k::T) where {T<:AbstractFloat}

Calculates the exact analytical static structure factor ``S(k)`` for a fluid interacting via
a Hard-Sphere + Two-Yukawa (SALR) potential using the Mean Spherical Approximation (MSA).

``V(r)/(k_B T) = -K_1 \\frac{e^{-z_1(r-1)}}{r} - K_2 \\frac{e^{-z_2(r-1)}}{r}, \\quad r > 1``

# Arguments
- `ϕ::T`: Particle volume fraction (``0 \\le \\phi < 1``).
- `K1::T`: Dimensionless energy parameter for the first Yukawa term. Positive for attraction, negative for repulsion.
- `z1::T`: Dimensionless inverse screening length for the first Yukawa term (``z_1 > 0``).
- `K2::T`: Dimensionless energy parameter for the second Yukawa term. Positive for attraction, negative for repulsion.
- `z2::T`: Dimensionless inverse screening length for the second Yukawa term (``z_2 > 0``).
- `k::T`: Dimensionless wavevector ``k = q\\sigma``.

# References
- Blum, L. (1980). Primitive conditions for the factorizable cluster expansion. J. Stat. Phys., 22(6), 661-672.
- Cummings, P. T., & Smith, W. R. (1979). Mol. Phys., 38(3), 997-1001.
- Liu, Y., Chen, W.-R., & Chen, S.-H. (2005). Cluster formation in two-Yukawa fluids. J. Chem. Phys., 122(4), 044507.
"""
function S_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k::T) where {T<:AbstractFloat}
    if ϕ < zero(T) || ϕ >= one(T)
        throw(DomainError(ϕ, "Volume fraction ϕ must be in [0, 1)."))
    end
    if z1 <= zero(T) || z2 <= zero(T)
        throw(DomainError((z1, z2), "Screening parameters z1, z2 must be strictly positive."))
    end
    if k < zero(T)
        throw(DomainError(k, "Wavevector k must be non-negative."))
    end

    if ϕ ≈ zero(T)
        return one(T)
    end

    # Hard-sphere limit
    if abs(K1) < eps(T) && abs(K2) < eps(T)
        return S_HS_PY(ϕ, k)
    end
    # Single-Yukawa limits
    if abs(K2) < eps(T)
        return S_Yukawa_MSA(ϕ, K1, z1, k)
    end
    if abs(K1) < eps(T)
        return S_Yukawa_MSA(ϕ, K2, z2, k)
    end

    c = _tysq_solve_physical_coefficients(ϕ, K1, z1, K2, z2)
    Ck_val = _tysq_Ck(c.a00, c.b00, c.v1, c.v2, k, ϕ, c.z1, c.K1, c.z2, c.K2)
    ρ = T(6.0) * ϕ / T(π)
    inv_S = one(T) - ρ * Ck_val
    return one(T) / inv_S
end

S_SALR_MSA(ϕ::Real, K1::Real, z1::Real, K2::Real, z2::Real, k::Real) =
    S_SALR_MSA(promote(Float64(ϕ), Float64(K1), Float64(z1), Float64(K2), Float64(z2), Float64(k))...)

function S_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T) where {T<:AbstractFloat}
    if abs(K1) < eps(T) && abs(K2) < eps(T)
        return k -> S_HS_PY(ϕ, Float64(k))
    elseif abs(K2) < eps(T)
        return S_Yukawa_MSA(ϕ, K1, z1)
    elseif abs(K1) < eps(T)
        return S_Yukawa_MSA(ϕ, K2, z2)
    end
    c = _tysq_solve_physical_coefficients(ϕ, K1, z1, K2, z2)
    return k -> begin
        Ck_val = _tysq_Ck(c.a00, c.b00, c.v1, c.v2, Float64(k), ϕ, c.z1, c.K1, c.z2, c.K2)
        ρ = 6.0 * ϕ / π
        1.0 / (1.0 - ρ * Ck_val)
    end
end
S_SALR_MSA(ϕ::Real, K1::Real, z1::Real, K2::Real, z2::Real) =
    S_SALR_MSA(promote(Float64(ϕ), Float64(K1), Float64(z1), Float64(K2), Float64(z2))...)

function S_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k_vec::AbstractVector{T}) where {T<:AbstractFloat}
    sk_func = S_SALR_MSA(ϕ, K1, z1, K2, z2)
    return [sk_func(k) for k in k_vec]
end

"""
    IS_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k::T) where {T<:AbstractFloat}

Calculates the inverse structure factor ``1 / S(k)`` for the Two-Yukawa (SALR) MSA model.
"""
function IS_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k::T) where {T<:AbstractFloat}
    return one(T) / S_SALR_MSA(ϕ, K1, z1, K2, z2, k)
end
IS_SALR_MSA(ϕ::Real, K1::Real, z1::Real, K2::Real, z2::Real, k::Real) =
    IS_SALR_MSA(promote(Float64(ϕ), Float64(K1), Float64(z1), Float64(K2), Float64(z2), Float64(k))...)

function IS_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T) where {T<:AbstractFloat}
    sk_func = S_SALR_MSA(ϕ, K1, z1, K2, z2)
    return k -> 1.0 / sk_func(k)
end
IS_SALR_MSA(ϕ::Real, K1::Real, z1::Real, K2::Real, z2::Real) =
    IS_SALR_MSA(promote(Float64(ϕ), Float64(K1), Float64(z1), Float64(K2), Float64(z2))...)

function IS_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k_vec::AbstractVector{T}) where {T<:AbstractFloat}
    isk_func = IS_SALR_MSA(ϕ, K1, z1, K2, z2)
    return [isk_func(k) for k in k_vec]
end

"""
    C_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k::T) where {T<:AbstractFloat}

Calculates the direct correlation function ``C(k) = \\frac{1 - 1/S(k)}{24 \\phi}`` for the Two-Yukawa (SALR) MSA model.
"""
function C_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k::T) where {T<:AbstractFloat}
    inv_S = IS_SALR_MSA(ϕ, K1, z1, K2, z2, k)
    return (one(T) - inv_S) / (T(24.0) * ϕ)
end
C_SALR_MSA(ϕ::Real, K1::Real, z1::Real, K2::Real, z2::Real, k::Real) =
    C_SALR_MSA(promote(Float64(ϕ), Float64(K1), Float64(z1), Float64(K2), Float64(z2), Float64(k))...)

function C_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T) where {T<:AbstractFloat}
    isk_func = IS_SALR_MSA(ϕ, K1, z1, K2, z2)
    return k -> (1.0 - isk_func(k)) / (24.0 * ϕ)
end
C_SALR_MSA(ϕ::Real, K1::Real, z1::Real, K2::Real, z2::Real) =
    C_SALR_MSA(promote(Float64(ϕ), Float64(K1), Float64(z1), Float64(K2), Float64(z2))...)

function C_SALR_MSA(ϕ::T, K1::T, z1::T, K2::T, z2::T, k_vec::AbstractVector{T}) where {T<:AbstractFloat}
    ck_func = C_SALR_MSA(ϕ, K1, z1, K2, z2)
    return [ck_func(k) for k in k_vec]
end

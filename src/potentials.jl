# --- Fourier Transformed Potentials ---

@doc """
    betaU_SW(temperature::T_num, lambda_range::T_num, k::T_num) where {T_num<:AbstractFloat}

Auxiliary function that computes the Fourier Transform of the attractive part of a
Square-Well (SW) potential, multiplied by β = 1/(k_B T_abs).
The result is β * u_attr_SW_tilde(k).

# Arguments
- `temperature::T_num`: Dimensionless temperature (often k_B * T_abs / ε, where ε is well depth).
- `lambda_range::T_num`: Range of the square well in terms of the hard-sphere diameter σ (i.e., λ*σ).
- `k::T_num`: Dimensionless wave vector (k = q*σ).

# Returns
- `T_num`: The value of β * u_attr_SW_tilde(k).
"""
function betaU_SW(temperature::T_num, lambda_range::T_num, k::T_num) where {T_num<:AbstractFloat}
    # This function calculates -β * u_SW_tilde(k) if u_SW is attractive (depth -ε)
    # and temperature is β = 1/ε (so 1/T_val is the well depth in units of k_B T_abs).
    # Or, if temperature includes ε, then it's -u_SW_tilde(k) / (k_B T_abs).
    # The sign convention depends on how the perturbation is defined.
    # The factor 4π comes from the spherical Fourier transform.

    k_abs = abs(k) # Use absolute value of k
    c_aux::T_num = zero(T_num)

    if k_abs > T_num(0.0750) # Threshold for using direct formula vs. Taylor expansion
        k2 = k_abs * k_abs
        λk = lambda_range * k_abs
        
        sink_val = sin(k_abs)
        cosk_val = cos(k_abs)
        sinλk_val = sin(λk)
        cosλk_val = cos(λk)
        
        # Term: (cos(k) - λ*cos(λk))/k³ + (λ*sin(λk) - sin(k))/(k²*λ)
        # The user's original form: c_aux = ((cosk - λ * cosλk) / k) + ((sinλk - sink) / k2); c_aux = c_aux / k
        # which is (cos(k) - λ*cos(λk))/k² + (λ*sin(λk) - sin(k))/k³
        # This form is common for the attractive part of the square-well potential.
        # Let's use the more standard form from literature if possible or stick to user's.
        # User's form: ( (cos(k) - λ*cos(λk))/k + (sin(λk) - sin(k))/k² ) / k
        # = (cos(k) - λ*cos(λk))/k² + (sin(λk) - sin(k))/k³
        term1_num = cosk_val - lambda_range * cosλk_val
        term2_num = sinλk_val - sink_val
        
        # Avoid division by zero if k_abs is extremely small but still > 0.0750 (unlikely but for safety)
        if k_abs < eps(T_num)
             c_aux = (one(T_num)/T_num(3.0)) * (lambda_range^3 - one(T_num)) # Fallback to k=0 limit
        else
            c_aux = (term1_num / k_abs + term2_num / k2) / k_abs
        end
    else 
        # Taylor expansion for small k
        λ3 = lambda_range^3
        λ5 = λ3 * lambda_range^2 # lambda_range^5
        k2_small = k_abs * k_abs
        
        c_aux = (one(T_num)/T_num(3.0)) * (λ3 - one(T_num)) - 
                (one(T_num)/T_num(30.0)) * (λ5 - one(T_num)) * k2_small
    end
    
    if temperature ≈ zero(T_num)
        # Avoid division by zero if temperature is zero.
        # Behavior depends on convention: if T is kBT/epsilon, then 1/T is large.
        # If T is just kBT and potential has depth epsilon, then this is -epsilon * c_aux / (kBT)
        # Assuming temperature is a positive value representing thermal energy scale.
        @warn "Temperature is zero in betaU_SW, result might be Inf or NaN."
        return sign(c_aux) * T_num(Inf) # Or handle as an error
    end
    
    return T_num(4.0) * T_num(π) * c_aux / temperature
end

@doc """
    betaU_Yukawa(amplitude::β_num, inv_screening_length::β_num, k::β_num) where {β_num<:AbstractFloat}

Auxiliary function to compute the Fourier Transform of a Yukawa potential,
multiplied by β = 1/(k_B T_abs).
The result is β * u_Yukawa_tilde(k).

The assumed dimensionless Yukawa potential in real space is u(r/σ) = (Amplitude) * exp(-z*(r/σ - 1)) / (r/σ) for r/σ >= 1.
The Fourier transform can take various forms depending on the precise definition
and normalization. The formula used here (from user) is:
βu_tilde(k) = Factor * Amplitude * (k*cos(k) + z*sin(k)) / (k*(k² + z²))
where k is qσ, z is inv_screening_length, and Factor is 4π.
The small k expansion is:
βu_tilde(k) ≈ Factor * Amplitude * ( (1/z) + 1 - (z²+9z+6)/(6z²) * k² )

Note: There's a discrepancy in the k -> 0 limit between the main formula and the small-k expansion
provided in the original code. The main formula for k->0 gives (1+z)/z², while the expansion gives (1+z)/z.
This implementation uses the formulas as provided.

# Arguments
- `amplitude::β_num`: Effective amplitude of the Yukawa potential (often includes β, e.g., βε).
- `inv_screening_length::β_num`: Dimensionless inverse screening length (z).
- `k::β_num`: Dimensionless wave vector (k = q*σ).

# Returns
- `β_num`: The value of β * u_Yukawa_tilde(k).

# References
[1] Yukawa, H. (1935). "On the interaction of elementary particles". Proc. Phys.-Math. Soc. Jpn. 17: 48.
    (Note: Specific forms for liquids often derive from this but might differ.)
"""
function betaU_Yukawa(amplitude::β_num, inv_screening_length::β_num, k::β_num) where {β_num<:AbstractFloat}
    k_abs = abs(k)
    z = inv_screening_length # Alias for clarity

    # Check for z=0 if it's problematic for the formulas
    if z ≈ zero(β_num) && k_abs ≈ zero(β_num) # Both zero, problematic for expansion
        @warn "Inverse screening length z and wavevector k are both near zero in betaU_Yukawa. Result may be Inf/NaN."
        # Return a value consistent with unscreened Coulomb if applicable, or handle as error.
        # For now, let it proceed to potential division by zero in expansion.
    elseif z ≈ zero(β_num) # Only z is zero (Coulomb limit)
         # The main formula becomes (k*cos(k))/(k*k^2) = cos(k)/k^2
         # The expansion for z=0 is ill-defined.
         if k_abs < β_num(0.075)
            @warn "Small-k expansion for betaU_Yukawa is ill-defined for z=0. Using main formula."
            # Fallback to main formula or a specific Coulomb FT if k_abs is also small
            if k_abs < eps(β_num)
                return β_num(Inf) # Or some other limit for Coulomb 1/k^2
            else
                return β_num(4.0) * β_num(π) * amplitude * (cos(k_abs) / k_abs^2)
            end
         else
            return β_num(4.0) * β_num(π) * amplitude * (cos(k_abs) / k_abs^2) # Simplified for z=0
         end
    end


    if k_abs < β_num(0.075)
        # Small k expansion (as provided by user)
        # Original: 4*π*A*((1/z)+1-(z*z+9*z+6)*k*k/(6*z*z))
        # Note: The k=0 limit of this is 4*π*A*(1/z + 1)
        term_const = one(β_num)/z + one(β_num)
        term_k2 = (z^2 + β_num(9.0)*z + β_num(6.0)) * k_abs^2 / (β_num(6.0)*z^2)
        return β_num(4.0) * β_num(π) * amplitude * (term_const - term_k2)
    else
        # Main formula
        # Original: 4*π*A*(k*cos(k) + z*sin(k))/(k*(k^2+z^2))
        # Note: The k=0 limit of this is 4*π*A*((1+z)/z^2)
        # There is a known discrepancy between the k=0 limits of the two branches.
        numerator = k_abs * cos(k_abs) + z * sin(k_abs)
        denominator = k_abs * (k_abs^2 + z^2)
        if denominator ≈ zero(β_num)
            @warn "Denominator is zero in betaU_Yukawa for k=$k_abs, z=$z. Result may be Inf/NaN."
            return sign(numerator) * β_num(Inf) # Or handle as error
        end
        return β_num(4.0) * β_num(π) * amplitude * numerator / denominator
    end
end


# Improved Lennard Jones kernel
function ILJ_kernel(ϵ, r0, β, m)
    n(r) = β + 4.0*(r/r0)^2
    return r -> ϵ * ( m/(n(r)-m)*(r0/r)^n(r) - n(r)/(n(r)-m)*(r0/r)^m )
end

# Specific induction kernel for the azide anion and carbon
function N3m_C_induction_kernel(α_C, Q_N_1, Q_N_2, Q_N_3)
    return (R_i1, R_i2, R_i3) -> -7200 * α_C * ( Q_N_1/(R_i1^2) + Q_N_2/(R_i2^2) + Q_N_3/(R_i3^2) )^2
end

# General electrostatic kernel
function electrostatic_kernel()
    return (Q1, Q2, r) -> 7200 * Q1 * Q2 / r
 end
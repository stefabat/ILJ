# Improved Lennard Jones kernel
# the energy unit of ϵ decideds the unit of this expression
# here ϵ is expected to be in meV
# out put is in meV units
function ILJ_kernel(ϵ, r0, β, m)
    n(r) = β + 4.0*(r/r0)^2
    return r -> ϵ * ( m/(n(r)-m)*(r0/r)^n(r) - n(r)/(n(r)-m)*(r0/r)^m )
end

# Specific induction kernel for the azide anion and carbon
# factor -1/2 in front is for the formula [Kaplan]
# α_C is supposed to be in Ang^3, charges in au units
# distances in Ang
# output is in meV units
function N3m_C_induction_kernel(α_C, Q_N_1, Q_N_2, Q_N_3, factor = 1.0)
    return (R_i1, R_i2, R_i3) -> -0.5*1000*toev/tobohr * α_C * factor * ( Q_N_1/(R_i1^2) + Q_N_2/(R_i2^2) + Q_N_3/(R_i3^2) )^2
end

# multiple dispatch
function N3m_C_induction_kernel(α_C, Q_N_1, Q_N_2, Q_N_3, β, r0)
    n(r) = β + 4.0*(r/r0)^2
    return (R_i1, R_i2, R_i3) -> -0.5*1000*toev/tobohr * α_C * n(R_i2)/(n(R_i2)-4) * ( Q_N_1/(R_i1^2) + Q_N_2/(R_i2^2) + Q_N_3/(R_i3^2) )^2
end

# General electrostatic kernel
# charges assumed to be in au units
# while distances assumed to be in Ang
# result is in meV units
function electrostatic_kernel()
    return (Q1, Q2, R_ij) -> (Q1 * Q2 / R_ij) * 1000*toev/tobohr
 end
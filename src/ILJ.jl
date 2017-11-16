# Improved Lennard Jones kernel
# the energy unit of ϵ decideds the unit of this expression
# here ϵ is expected to be in meV
# out put is in meV units
function ILJ_kernel(ϵ, r0, β, m)
    if ϵ == 0
        return r -> 0.0
    end
    n(r) = β + 4.0*(r/r0)^2
    return r -> ϵ * ( m/(n(r)-m)*(r0/r)^n(r) - n(r)/(n(r)-m)*(r0/r)^m )
end

# Specific induction kernel for the azide anion and carbon
# factor -1/2 in front is for the formula [Kaplan]
# α is supposed to be in Ang^3, charges in au units
# distances in Ang
# output is in meV units
function N3m_induction_kernel(α, Q_N_1, Q_N_2, Q_N_3, factor = 1.0)
    return (R_i1, R_i2, R_i3) -> -0.5*1000*toev/tobohr * α * factor * ( Q_N_1/(R_i1^2) + Q_N_2/(R_i2^2) + Q_N_3/(R_i3^2) )^2
end

# multiple dispatch
function N3m_induction_kernel(α, Q_N_1, Q_N_2, Q_N_3, β, r0)
    n(r) = β + 4.0*(r/r0)^2
    return (R_i1, R_i2, R_i3) -> -0.5*1000*toev/tobohr * α * n(R_i2)/(n(R_i2)-4) * ( Q_N_1/(R_i1^2) + Q_N_2/(R_i2^2) + Q_N_3/(R_i3^2) )^2
end

# General electrostatic kernel
# charges assumed to be in au units
# while distances assumed to be in Ang
# result is in meV units
function electrostatic_kernel()
    return (Q1, Q2, R_ij) -> (Q1 * Q2 / R_ij) * 1000*toev/tobohr
 end

function Rm(α1_eff,α2_eff)
    return 1.767 * (α1_eff^(1/3) + α2_eff^(1/3))/(α1_eff*α2_eff)^0.095
end

function C6eff(α1,N1_eff,α2,N2_eff)
    return 15700 * α1 * α2 / (sqrt(α1/N1_eff) + sqrt(α2/N2_eff))
end

function epsilon(R,C6)
    return 0.720 * C6/R^6
end


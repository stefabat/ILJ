# Improved Lennard Jones kernel

function ILJ_kernel(ϵ, r0, β, m)
    n(r) = β + 4.0*(r/r0)^2
    return r -> ϵ * [ m/(n(r)-m)*(r0/r)^n(r) - n(r)/(n(r)-m)*(r0/r)^m ]
end

function N3_induction_kernel(r0, β, q1, q2, q3, α_c)
    n(r) = β + 4.0*(r/r0)^2
    return (R1,R2,R3) -> -7200 * n(R2)/(n(R2) - 4) * ( (q1/(R1)^2) + (q2/(R2)^2) + (q3/(R3)^2))^2 * α_c
end

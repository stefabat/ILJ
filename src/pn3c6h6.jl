
include("ILJ.jl")
include("xyz_parser.jl")

# in kcal/mol
ϵ_C_N = 0.8443; rm_C_N = 4.220;      # Carbon to external nitrogens
ϵ_H_N = 0.4243; rm_H_N = 3.730;
β = 8.0
m = 4.0
q_Next = -0.56; q_Nint = +0.12; # in au units
α_C = 1.20 # in Ang^3
α_H = 0.38
    
# Read the two fragments
frag_1 = readxyz("data/geometries/ch6.xyz")
frag_2 = readxyz("data/geometries/n3.xyz")

# Assign charges to fragments
frag_1.charges = repeat([-0.1;0.1], inner=6)
frag_2.charges = [q_Next;q_Nint;q_Next]

# ILJ potentials
ILJ_C_N = ILJ_kernel(ϵ_C_N, rm_C_N, β, m)
ILJ_H_N = ILJ_kernel(ϵ_H_N, rm_H_N, β, m)

# Electrostatic potential
els = electrostatic_kernel()

# Initialize potential sums
V_ILJ_C_N  = 0.0;
V_ILJ_H_N  = 0.0;
V_els_C_N1 = 0.0; V_els_C_N2 = 0.0; V_els_C_N3 = 0.0;
V_els_H_N1 = 0.0; V_els_H_N2 = 0.0; V_els_H_N3 = 0.0;

# Conversion factor meV to kcal/mol
mevtokcal = 0.023060554446846782

# Loop over fragment 1
for i in range(1,frag_1.n_atoms)
    ri1 = vecnorm(frag_1.coords[i,:] - frag_2.coords[1,:])
    ri2 = vecnorm(frag_1.coords[i,:] - frag_2.coords[2,:])
    ri3 = vecnorm(frag_1.coords[i,:] - frag_2.coords[3,:])
    if frag_1.types[i] == "C"
        V_els_C_N1 += els(frag_1.charges[i], frag_2.charges[1], ri1)
        V_els_C_N2 += els(frag_1.charges[i], frag_2.charges[2], ri2)
        V_els_C_N3 += els(frag_1.charges[i], frag_2.charges[3], ri3)
        V_ILJ_C_N  += ILJ_C_N(ri2)
    elseif frag_1.types[i] == "H"
        V_els_H_N1 += els(frag_1.charges[i], frag_2.charges[1], ri1)
        V_els_H_N2 += els(frag_1.charges[i], frag_2.charges[2], ri2)
        V_els_H_N3 += els(frag_1.charges[i], frag_2.charges[3], ri3)
        V_ILJ_H_N  += ILJ_H_N(ri2)
    end
end

# in kcal
V_ILJ_tot = V_ILJ_C_N + V_ILJ_H_N
V_els_tot = (V_els_C_N1 + V_els_C_N2 + V_els_C_N3 + V_els_H_N1 + V_els_H_N2 + V_els_H_N3)*tokcal
V_int_tot = V_ILJ_tot  + V_els_tot

##### PRINTING OUTPUT #####
println("\n\tParameters:     β = ",β,"    m = ",m,"    α_C = ",α_C,"    α_H = ",α_H)
println("\t\t\tr0_C_N = ",round(rm_C_N,3),"\t     ϵ_C_N = ",round(ϵ_C_N,3))
println("\t\t\tr0_H_N = ",round(rm_H_N,3),"\t     ϵ_H_N = ",round(ϵ_H_N,3))

@printf("\n%46s %13.5f [meV]\n","Total Improved Lennard-Jones term:",V_ILJ_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total Improved Lennard-Jones term:",V_ILJ_tot)

@printf("%46s %13.5f [meV]\n","Total electrostatic term:",V_els_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total electrostatic term:",V_els_tot)

@printf("%46s %13.5f [meV]\n","Total interaction energy:",V_int_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total interaction energy:",V_int_tot)


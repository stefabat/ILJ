
include("ILJ.jl")
include("xyz_parser.jl")

# in meV
#ϵ_C_Next = 5.205; rm_C_Next = 3.994;      # Carbon to external nitrogens
#ϵ_C_Nint = 3.536; rm_C_Nint = 3.818;      # Carbon to internal nitrogen
#ϵ_H_Next = 2.827; rm_H_Next = 3.644;
#ϵ_H_Nint = 2.431; rm_H_Nint = 3.348;
# in kcal/mol
ϵ_C_Next = 0.1210; rm_C_Next = 3.994;      # Carbon to external nitrogens
ϵ_C_Nint = 0.0822; rm_C_Nint = 3.818;      # Carbon to internal nitrogen
ϵ_H_Next = 0.0657; rm_H_Next = 3.644;
ϵ_H_Nint = 0.0565; rm_H_Nint = 3.348;
β = 8.0
m = 6.0
q_Next = -0.56; q_Nint = +0.12; # in au units
α_C = 1.136 # in Ang^3
α_H = 0.38
    
# Read the two fragments
frag_1 = readxyz("data/geometries/ch6.xyz")
frag_2 = readxyz("data/geometries/n3.xyz")

# Assign charges to fragments
frag_1.charges = repeat([-0.1;0.1], inner=6)
frag_2.charges = [q_Next;q_Nint;q_Next]

# ILJ potentials
ILJ_C_Next = ILJ_kernel(ϵ_C_Next, rm_C_Next, β, m)
ILJ_C_Nint = ILJ_kernel(ϵ_C_Nint, rm_C_Nint, β, m)
ILJ_H_Next = ILJ_kernel(ϵ_H_Next, rm_H_Next, β, m)
ILJ_H_Nint = ILJ_kernel(ϵ_H_Nint, rm_H_Nint, β, m)

# Induction potential specific for azide-carbon interaction
ind_C_N3m = N3m_induction_kernel(α_C, q_Next, q_Nint, q_Next)
ind_H_N3m = N3m_induction_kernel(α_H, q_Next, q_Nint, q_Next)

# Electrostatic potential
els = electrostatic_kernel()

# Initialize potential sums
V_ILJ_C_N1 = 0.0; V_ILJ_C_N2 = 0.0; V_ILJ_C_N3 = 0.0;
V_ILJ_H_N1 = 0.0; V_ILJ_H_N2 = 0.0; V_ILJ_H_N3 = 0.0;
V_ind_C_N  = 0.0; V_ind_H_N  = 0.0;
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
        V_ILJ_C_N1 += ILJ_C_Next(ri1)
        V_ILJ_C_N2 += ILJ_C_Nint(ri2)
        V_ILJ_C_N3 += ILJ_C_Next(ri3)
        V_ind_C_N  += ind_C_N3m(ri1, ri2, ri3)
    elseif frag_1.types[i] == "H"
        V_els_H_N1 += els(frag_1.charges[i], frag_2.charges[1], ri1)
        V_els_H_N2 += els(frag_1.charges[i], frag_2.charges[2], ri2)
        V_els_H_N3 += els(frag_1.charges[i], frag_2.charges[3], ri3)
        V_ILJ_H_N1 += ILJ_H_Next(ri1)
        V_ILJ_H_N2 += ILJ_H_Nint(ri2)
        V_ILJ_H_N3 += ILJ_H_Next(ri3)
        V_ind_H_N  += ind_H_N3m(ri1, ri2, ri3)
    end
end

# in kcal
V_ILJ_tot = V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3 + V_ILJ_H_N1 + V_ILJ_H_N2 + V_ILJ_H_N3
V_els_tot = (V_els_C_N1 + V_els_C_N2 + V_els_C_N3 + V_els_H_N1 + V_els_H_N2 + V_els_H_N3)*tokcal
V_ind_tot = (V_ind_C_N  + V_ind_H_N)*tokcal
V_int_tot = V_ILJ_tot  + V_els_tot  + V_ind_tot

##### PRINTING OUTPUT #####
println("\n\tParameters:     β = ",β,"    m = ",m,"    α_C = ",α_C,"    α_H = ",α_H)
println("\t\t\tr0_C_Ni = ",round(rm_C_Nint,3),"\t     ϵ_C_Ni = ",round(ϵ_C_Nint,3))
println("\t\t\tr0_C_Ne = ",round(rm_C_Next,3),"\t     ϵ_C_Ne = ",round(ϵ_C_Next,3))
println("\t\t\tr0_H_Ni = ",round(rm_H_Nint,3),"\t     ϵ_H_Ni = ",round(ϵ_H_Nint,3))
println("\t\t\tr0_H_Ne = ",round(rm_H_Next,3),"\t     ϵ_H_Ne = ",round(ϵ_H_Next,3))

@printf("\n%46s %13.5f [meV]\n","Total Improved Lennard-Jones term:",V_ILJ_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total Improved Lennard-Jones term:",V_ILJ_tot)

@printf("%46s %13.5f [meV]\n","Total electrostatic term:",V_els_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total electrostatic term:",V_els_tot)

@printf("%46s %13.5f [meV]\n","Total induction term:",V_ind_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total induction term:",V_ind_tot)

@printf("%46s %13.5f [meV]\n","Total interaction energy:",V_int_tot/mevtokcal)
@printf("%46s %13.5f [kcal/mol]\n\n","Total interaction energy:",V_int_tot)


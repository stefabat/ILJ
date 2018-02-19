
include("ILJ.jl")
include("xyz_parser.jl")

ϵ_C_N_ext = 5.205; r0_C_N_ext = 3.994;      # Carbon to external nitrogens
ϵ_C_N_int = 3.536; r0_C_N_int = 3.818;      # Carbon to internal nitrogen
ϵ_H_N_ext = 2.827; r0_H_N_ext = 3.644;
ϵ_H_N_int = 2.431; r0_H_N_int = 3.348;
β = 8.0
m = 6.0
Q_N_ext = -0.56; Q_N_int = +0.12; # in au units
α_C = 1.136 # in Ang^3
α_H = 0.38

# Read the two fragments
frag_1 = readxyz("data/geometries/ch6.xyz")
frag_2 = readxyz("data/geometries/n3.xyz")

# ILJ potentials
ILJ_C_N_ext = ILJ_kernel(ϵ_C_N_ext, r0_C_N_ext, β, m)
ILJ_C_N_int = ILJ_kernel(ϵ_C_N_int, r0_C_N_int, β, m)
ILJ_H_N_ext = ILJ_kernel(ϵ_H_N_ext, r0_H_N_ext, β, m)
ILJ_H_N_int = ILJ_kernel(ϵ_H_N_int, r0_H_N_int, β, m)

# Induction potential specific for azide-carbon interaction
ind_C_N3m = N3m_induction_kernel(α_C, Q_N_ext, Q_N_int, Q_N_ext, 1.0)
ind_H_N3m = N3m_induction_kernel(α_H, Q_N_ext, Q_N_int, Q_N_ext, 1.0)

# Initialize potential sums
V_ILJ_C_N1 = 0.0; V_ILJ_C_N2 = 0.0; V_ILJ_C_N3 = 0.0;
V_ILJ_H_N1 = 0.0; V_ILJ_H_N2 = 0.0; V_ILJ_H_N3 = 0.0;
V_ind_C_N  = 0.0; V_ind_H_N  = 0.0;

# Conversion factor meV to kcal/mol
mevtokcal = 0.023060554446846782

# Loop over fragment 1
for i in range(1,frag_1.n_atoms)
    R1 = vecnorm(frag_1.coords[i,:] - frag_2.coords[1,:])
    R2 = vecnorm(frag_1.coords[i,:] - frag_2.coords[2,:])
    R3 = vecnorm(frag_1.coords[i,:] - frag_2.coords[3,:])
    if frag_1.types[i] == "C"
        V_ILJ_C_N1 += ILJ_C_N_ext(R1)
        V_ILJ_C_N2 += ILJ_C_N_int(R2)
        V_ILJ_C_N3 += ILJ_C_N_ext(R3)
        V_ind_C_N  += ind_C_N3m(R1, R2, R3)
    elseif frag_1.types[i] == "H"
        V_ILJ_H_N1 += ILJ_H_N_ext(R1)
        V_ILJ_H_N2 += ILJ_H_N_int(R2)
        V_ILJ_H_N3 += ILJ_H_N_ext(R3)
        V_ind_H_N  += ind_H_N3m(R1, R2, R3)
    end
end

# Sum everything in meV
V_ILJ_tot = V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3 + V_ILJ_H_N1 + V_ILJ_H_N2 + V_ILJ_H_N3
V_ind_tot = V_ind_C_N  + V_ind_H_N
V_int_tot = V_ILJ_tot  + V_ind_tot

# Transform in kcal/mol
V_ILJ_tot_kcal = V_ILJ_tot * mevtokcal
V_ind_tot_kcal = V_ind_tot * mevtokcal
V_int_tot_kcal = V_int_tot * mevtokcal

##### PRINTING OUTPUT #####
println("\n\tParameters:     β = ",β,"    m = ",m,"    α_C = ",α_C,"    α_H = ",α_H)
println("\t\t\tr0_C_Ni = ",round(r0_C_N_int,3),"\t     ϵ_C_Ni = ",round(ϵ_C_N_int,3))
println("\t\t\tr0_C_Ne = ",round(r0_C_N_ext,3),"\t     ϵ_C_Ne = ",round(ϵ_C_N_ext,3))
println("\t\t\tr0_H_Ni = ",round(r0_H_N_int,3),"\t     ϵ_H_Ni = ",round(ϵ_H_N_int,3))
println("\t\t\tr0_H_Ne = ",round(r0_H_N_ext,3),"\t     ϵ_H_Ne = ",round(ϵ_H_N_ext,3))

@printf("\n%46s %8.2f [meV]\n","Total Improved Lennard-Jones term:",V_ILJ_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total Improved Lennard-Jones term:",V_ILJ_tot_kcal)

@printf("%46s %8.2f [meV]\n","Total induction term:",V_ind_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total induction term:",V_ind_tot_kcal)

@printf("%46s %8.2f [meV]\n","Total interaction energy:",V_int_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total interaction energy:",V_int_tot_kcal)


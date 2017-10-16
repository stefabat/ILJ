# Simple script to compute ILJ intearction energies
include("ILJ.jl")
include("xyz_parser.jl")

# Define ILJ parameters, ϵ is in meV and r is in Ang
ϵ_C_N_ext = 5.205; r0_C_N_ext = 3.994;      # Carbon to external nitrogens
ϵ_C_N_int = 3.536; r0_C_N_int = 3.818;      # Carbon to internal nitrogen
ϵ_H_N_ext = 2.827; r0_H_N_ext = 3.644;
ϵ_H_N_int = 2.431; r0_H_N_int = 3.348;
β = 8.0
m = 6.0
Q_N_ext = -0.56; Q_N_int = +0.12; # in au units
α_C = 1.2 # in Ang^3
α_H = 0.38
#α_C = 1.334 # in Ang^3
#α_H = 0.496

# Read the two fragments
frag_1 = readxyz("data/geometries/t5511_b97d3_cc-pvtz.xyz")
frag_2 = readxyz("data/geometries/n3_b97d3_aug-cc-pvtz.xyz")

# Parse atomic charges of fragment_1 and set them
frag_1.charges = readdlm("data/charges/t5511_b97d3_cc-pvtz_npa.dat")

# Set partial atomic charges to azide anion
frag_2.charges = [Q_N_ext Q_N_int Q_N_ext]

# ILJ potentials
ILJ_C_N_ext = ILJ_kernel(ϵ_C_N_ext, r0_C_N_ext, β, m)
ILJ_C_N_int = ILJ_kernel(ϵ_C_N_int, r0_C_N_int, β, m)
ILJ_H_N_ext = ILJ_kernel(ϵ_H_N_ext, r0_H_N_ext, β, m)
ILJ_H_N_int = ILJ_kernel(ϵ_H_N_int, r0_H_N_int, β, m)

# Induction potential specific for azide-carbon interaction
ind_C_N3m = N3m_induction_kernel(α_C, Q_N_ext, Q_N_int, Q_N_ext, 1.0)
ind_H_N3m = N3m_induction_kernel(α_H, Q_N_ext, Q_N_int, Q_N_ext, 1.0)
# ind_C_N3m = N3m_C_induction_kernel(α_C, Q_N_ext, Q_N_int, Q_N_ext, β_C_N, r0_C_N_int)

# Electrostatic potential
els = electrostatic_kernel()

# Initialize potential sums
V_ILJ_C_N1 = 0.0; V_ILJ_C_N2 = 0.0; V_ILJ_C_N3 = 0.0;
V_ILJ_H_N1 = 0.0; V_ILJ_H_N2 = 0.0; V_ILJ_H_N3 = 0.0;
V_ind_C_N  = 0.0; V_ind_H_N  = 0.0;
V_els_C_N1 = 0.0; V_els_C_N2 = 0.0; V_els_C_N3 = 0.0;
V_els_H_N1 = 0.0; V_els_H_N2 = 0.0; V_els_H_N3 = 0.0;

# Loop over fragment 1
for i in range(1,frag_1.n_atoms)
    R1 = vecnorm(frag_1.coords[i,:] - frag_2.coords[1,:])
    R2 = vecnorm(frag_1.coords[i,:] - frag_2.coords[2,:])
    R3 = vecnorm(frag_1.coords[i,:] - frag_2.coords[3,:])
    if frag_1.types[i] == "C"
        V_els_C_N1 += els(frag_1.charges[i], frag_2.charges[1], R1)
        V_els_C_N2 += els(frag_1.charges[i], frag_2.charges[2], R2)
        V_els_C_N3 += els(frag_1.charges[i], frag_2.charges[3], R3)
        V_ILJ_C_N1 += ILJ_C_N_ext(R1)
        V_ILJ_C_N2 += ILJ_C_N_int(R2)
        V_ILJ_C_N3 += ILJ_C_N_ext(R3)
        V_ind_C_N  += ind_C_N3m(R1, R2, R3)
    elseif frag_1.types[i] == "H"
        V_els_H_N1 += els(frag_1.charges[i], frag_2.charges[1], R1)
        V_els_H_N2 += els(frag_1.charges[i], frag_2.charges[2], R2)
        V_els_H_N3 += els(frag_1.charges[i], frag_2.charges[3], R3)
        V_ILJ_H_N1 += ILJ_H_N_ext(R1)
        V_ILJ_H_N2 += ILJ_H_N_int(R2)
        V_ILJ_H_N3 += ILJ_H_N_ext(R3)
        V_ind_H_N  += ind_H_N3m(R1, R2, R3)
    end
end

# Conversion factor meV to kcal/mol
mevtokcal = 0.023060554446846782

# Sum everything in meV
V_ILJ_tot = V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3 + V_ILJ_H_N1 + V_ILJ_H_N2 + V_ILJ_H_N3
V_els_tot = V_els_C_N1 + V_els_C_N2 + V_els_C_N3 + V_els_H_N1 + V_els_H_N2 + V_els_H_N3
V_ind_tot = V_ind_C_N  + V_ind_H_N
V_int_tot = V_ILJ_tot  + V_els_tot  + V_ind_tot

# Transform in kcal/mol
V_ILJ_tot_kcal = V_ILJ_tot * mevtokcal
V_els_tot_kcal = V_els_tot * mevtokcal
V_ind_tot_kcal = V_ind_tot * mevtokcal
V_int_tot_kcal = V_int_tot * mevtokcal

##### PRINTING OUTPUT #####
println("\n\tParameters:    β = ",β,"    m = ",m,"    α_C = ",α_C,"    α_H = ",α_H)

@printf("\n%46s %8.2f [meV]\n","Total Improved Lennard-Jones term:",V_ILJ_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total Improved Lennard-Jones term:",V_ILJ_tot_kcal)

@printf("%46s %8.2f [meV]\n","Total electrostatic term:",V_els_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total electrostatic term:",V_els_tot_kcal)

@printf("%46s %8.2f [meV]\n","Total induction term:",V_ind_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total induction term:",V_ind_tot_kcal)

@printf("%46s %8.2f [meV]\n","Total interaction energy:",V_int_tot)
@printf("%46s %8.2f [kcal/mol]\n\n","Total interaction energy:",V_int_tot_kcal)

# DEBUG printing
# println("\n----------------------------------\nDEBUG printing [in kcal/mol]\n----------------------------------")
# println("ILJ C N1: ",round(V_ILJ_C_N1*mevtokcal,2))
# println("ILJ C N2: ",round(V_ILJ_C_N2*mevtokcal,2))
# println("ILJ C N3: ",round(V_ILJ_C_N3*mevtokcal,2))
# println("")
# println("ILJ H N1: ",round(V_ILJ_H_N1*mevtokcal,2))
# println("ILJ H N2: ",round(V_ILJ_H_N2*mevtokcal,2))
# println("ILJ H N3: ",round(V_ILJ_H_N3*mevtokcal,2))
# println("\n")
# println("els C N1: ",round(V_els_C_N1*mevtokcal,2))
# println("els C N2: ",round(V_els_C_N2*mevtokcal,2))
# println("els C N3: ",round(V_els_C_N3*mevtokcal,2))
# println("")
# println("els H N1: ",round(V_els_H_N1*mevtokcal,2))
# println("els H N2: ",round(V_els_H_N2*mevtokcal,2))
# println("els H N3: ",round(V_els_H_N3*mevtokcal,2))
# println("\n")
# println("ind C Nx: ",round(V_ind_C_N *mevtokcal,2))
# println("ind H Nx: ",round(V_ind_H_N *mevtokcal,2))

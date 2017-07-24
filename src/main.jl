# Simple script to compute ILJ intearction energies
include("ILJ.jl")
include("xyz_parser.jl")

# Define ILJ parameters
ϵ_C_N_ext = 5.205; r0_C_N_ext = 3.994;      # Carbon to external nitrogens
ϵ_C_N_int = 3.536; r0_C_N_int = 3.818;      # Carbon to internal nitrogen
β_C_N = 8.0
m_C_N = 6.0
Q_N_ext = -0.56; Q_N_int = +0.12;
α_C = 1.2

# Read the two fragments
frag_1 = readxyz("data/t553_b97d3_cc-pvtz.xyz")
frag_2 = readxyz("data/n3_b97d3_aug-cc-pvtz.xyz")

# Set partial atomic charges to azide anion
frag_2.charges = [Q_N_ext Q_N_int Q_N_ext]

# Parse atomic charges of fragment_1 and set them
frag_1.charges = atomic_charges("data/t553_b97d3_cc-pvdz.log")
frag_1.charges[100] = 0.22998

# ILJ potentials
ILJ_C_N_ext = ILJ_kernel(ϵ_C_N_ext, r0_C_N_ext, β_C_N, m_C_N)
ILJ_C_N_int = ILJ_kernel(ϵ_C_N_int, r0_C_N_int, β_C_N, m_C_N)

# Induction potential specific for azide-carbon interaction
ind_C_N3m = N3m_C_induction_kernel(α_C, Q_N_ext, Q_N_int, Q_N_ext)

# Electrostatic potential
els = electrostatic_kernel()

# Initialize potential sums
V_ILJ_C_N1 = 0.0; V_ILJ_C_N2 = 0.0; V_ILJ_C_N3 = 0.0
V_ind_C_N  = 0.0
V_els_C_N1 = 0.0; V_els_C_N2 = 0.0; V_els_C_N3 = 0.0
V_els_H_N1 = 0.0; V_els_H_N2 = 0.0; V_els_H_N3 = 0.0

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
    end
end

# DEBUG printing
println("ILJ C N1: ",V_ILJ_C_N1)
println("ILJ C N2: ",V_ILJ_C_N2)
println("ILJ C N3: ",V_ILJ_C_N3)
println("ind C Nx: ",V_ind_C_N)
println("els C N1: ",V_els_C_N1)
println("els C N2: ",V_els_C_N2)
println("els C N3: ",V_els_C_N3)
println("els H N1: ",V_els_H_N1)
println("els H N2: ",V_els_H_N2)
println("els H N3: ",V_els_H_N3)

# Sum everything in meV
V_ILJ_tot = (V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3)
V_els_C_tot = (V_els_C_N1 + V_els_C_N2 + V_els_C_N3)
V_els_H_tot = (V_els_H_N1 + V_els_H_N2 + V_els_H_N3)
V_int_tot = V_ILJ_tot + V_els_C_tot + V_els_H_tot + V_ind_C_N

# Transform in kcal/mol
V_ILJ_tot_kcal = V_ILJ_tot*0.02306
V_els_C_tot_kcal = V_els_C_tot*0.02306
V_els_H_tot_kcal = V_els_H_tot*0.02306
V_ind_C_N_kcal = V_ind_C_N*0.02306
V_int_tot_kcal = V_int_tot*0.02306

#println("Total Improved Lennard-Jones interaction energy: ",round(V_ILJ_tot_kcal,2)," [kcal/mol]")
@printf("\n%50s %8.2f [meV]\n","Total Improved Lennard-Jones interaction energy:",V_ILJ_tot)
@printf("%50s %8.2f [kcal/mol]\n\n","Total Improved Lennard-Jones interaction energy:",V_ILJ_tot_kcal)

@printf("\n%50s %8.2f [meV]\n","Total C-N electrostatic interaction energy:",V_els_C_tot)
@printf("%50s %8.2f [kcal/mol]\n\n","Total C-N electrostatic interaction energy:",V_els_C_tot_kcal)

@printf("\n%50s %8.2f [meV]\n","Total H-N electrostatic interaction energy:",V_els_H_tot)
@printf("%50s %8.2f [kcal/mol]\n\n","Total H-N electrostatic interaction energy:",V_els_H_tot_kcal)

#println("Induction term interaction energy: ",round(V_ind_C_N_kcal,2)," [kcal/mol]")
@printf("%50s %8.2f [meV]\n","Induction term interaction energy:",V_ind_C_N)
@printf("%50s %8.2f [kcal/mol]\n\n","Induction term interaction energy:",V_ind_C_N_kcal)

#println("Total interaction energy: ",round(V_int_tot,2)," [kcal/mol]")
@printf("%50s %8.2f [meV]\n","Total interaction energy:",V_int_tot)
@printf("%50s %8.2f [kcal/mol]\n\n","Total interaction energy:",V_int_tot_kcal)


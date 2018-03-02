# Simple script to compute ILJ intearction energies
include("ILJ.jl")
include("xyz_parser.jl")

# Compute ILJ parameters
#α_C = 1.200; α_C_eff = 2.50; N_C_eff = 3.3
α_C = 1.136; α_C_eff = 2.414; N_C_eff = 3.3
α_H = 0.38;  α_H_eff = 0.40;  N_H_eff = 1.0
α_Next = 1.9 ; α_Next_eff = 2.20; N_Next_eff = 5.0
α_Nint = 0.85; α_Nint_eff = 1.00; N_Nint_eff = 4.2

# in Ang
rm_C_Next = Rm(α_C_eff,α_Next_eff)
rm_C_Nint = Rm(α_C_eff,α_Nint_eff)
rm_H_Next = Rm(α_H_eff,α_Next_eff)
rm_H_Nint = Rm(α_H_eff,α_Nint_eff)

C6_C_Next = C6eff(α_C,N_C_eff,α_Next,N_Next_eff)
C6_C_Nint = C6eff(α_C,N_C_eff,α_Nint,N_Nint_eff)
C6_H_Next = C6eff(α_H,N_H_eff,α_Next,N_Next_eff)
C6_H_Nint = C6eff(α_H,N_H_eff,α_Nint,N_Nint_eff)

# in meV
ϵ_C_Next = epsilon(rm_C_Next,C6_C_Next) 
ϵ_C_Nint = epsilon(rm_C_Nint,C6_C_Nint) 
ϵ_H_Next = epsilon(rm_H_Next,C6_H_Next) 
ϵ_H_Nint = epsilon(rm_H_Nint,C6_H_Nint) 

# Define ILJ parameters, ϵ is in meV and r is in Ang
#ϵ_C_N_ext = 5.205; r0_C_N_ext = 3.994;      # Carbon to external nitrogens
#ϵ_C_N_int = 3.536; r0_C_N_int = 3.818;      # Carbon to internal nitrogen
#ϵ_H_N_ext = 2.827; r0_H_N_ext = 3.644;
#ϵ_H_N_int = 2.431; r0_H_N_int = 3.348;
β = 8.0
m = 6.0
q_Next = -0.56; q_Nint = +0.12; # in au units
#α_C = 1.2 # in Ang^3
#α_H = 0.38

# Read the two fragments
frag_1 = readxyz("data/geometries/t553_b97d3_cc-pvtz.xyz")
frag_2 = readxyz("data/geometries/n3_b97d3_aug-cc-pvtz.xyz")

# Parse atomic charges of fragment_1 and set them
frag_1.charges = readdlm("data/charges/cc-pvtz/t553_b97d3_cc-pvtz_npa.dat")

# Set partial atomic charges to azide anion
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

# Loop over fragment 1
for i in range(1,frag_1.n_atoms)
    ri1 = vecnorm(frag_1.coords[i,:] - frag_2.coords[1,:])
    ri2 = vecnorm(frag_1.coords[i,:] - frag_2.coords[2,:])
    ri3 = vecnorm(frag_1.coords[i,:] - frag_2.coords[3,:])
    println(ri1,"  ",ri2,"  ",ri3)
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

# Conversion factor meV to kcal/mol
mevtokcal = 0.023060554446846782

# Sum everything in meV
V_ILJ_tot = V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3 + V_ILJ_H_N1 + V_ILJ_H_N2 + V_ILJ_H_N3
V_els_tot = (V_els_C_N1 + V_els_C_N2 + V_els_C_N3 + V_els_H_N1 + V_els_H_N2 + V_els_H_N3)*tomev
V_ind_tot = (V_ind_C_N  + V_ind_H_N)*tomev
V_int_tot = V_ILJ_tot  + V_els_tot  + V_ind_tot

# Transform in kcal/mol
V_ILJ_tot_kcal = V_ILJ_tot * mevtokcal
V_els_tot_kcal = V_els_tot * mevtokcal
V_ind_tot_kcal = V_ind_tot * mevtokcal
V_int_tot_kcal = V_int_tot * mevtokcal

##### PRINTING OUTPUT #####
println("\n\tParameters:     β = ",β,"    m = ",m,"    α_C = ",α_C,"    α_H = ",α_H)
println("\t\t\tr0_C_Ni = ",round(rm_C_Nint,3),"\t     ϵ_C_Ni = ",round(ϵ_C_Nint,3))
println("\t\t\tr0_C_Ne = ",round(rm_C_Next,3),"\t     ϵ_C_Ne = ",round(ϵ_C_Next,3))
println("\t\t\tr0_H_Ni = ",round(rm_H_Nint,3),"\t     ϵ_H_Ni = ",round(ϵ_H_Nint,3))
println("\t\t\tr0_H_Ne = ",round(rm_H_Next,3),"\t     ϵ_H_Ne = ",round(ϵ_H_Next,3))

@printf("\n%46s %11.5f [meV]\n","Total Improved Lennard-Jones term:",V_ILJ_tot)
@printf("%46s %11.5f [kcal/mol]\n\n","Total Improved Lennard-Jones term:",V_ILJ_tot_kcal)

@printf("%46s %11.5f [meV]\n","Total electrostatic term:",V_els_tot)
@printf("%46s %11.5f [kcal/mol]\n\n","Total electrostatic term:",V_els_tot_kcal)

@printf("%46s %11.5f [meV]\n","Total induction term:",V_ind_tot)
@printf("%46s %11.5f [kcal/mol]\n\n","Total induction term:",V_ind_tot_kcal)

@printf("%46s %11.5f [meV]\n","Total interaction energy:",V_int_tot)
@printf("%46s %11.5f [kcal/mol]\n\n","Total interaction energy:",V_int_tot_kcal)

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

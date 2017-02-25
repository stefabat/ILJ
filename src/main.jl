#script to execute int energy
include("ILJ.jl")
include("xyz_parser.jl")

# Define parameters
ϵ_C_N1 = 5.205; ϵ_C_N2 = 3.536;
r0_C_N1 = 3.994; r0_C_N2 = 3.818;
β_C_N1 = 9.0; β_C_N2 = 9.0;
m_C_N1 = 6.0; m_C_N2 = 6.0;
q1 = -0.56; q2 = +0.12; q3 = -0.56;
α_c = 1.2

# Define ILJ potentials
ILJ_C_N1 = ILJ_kernel(ϵ_C_N1, r0_C_N1, β_C_N1, m_C_N1)
ILJ_C_N2 = ILJ_kernel(ϵ_C_N2, r0_C_N2, β_C_N2, m_C_N2)

# Define induction potential
ind_C_N = N3_induction_kernel(r0_C_N2, β_C_N2,q1, q2, q3, α_c)

# Read geometry
atom_types,atom_coords = xyz_parser("data/n3_t555_1.xyz")

C_idx  = find(x -> x == "C",atom_types)
N1_idx = find(x -> x == "N1",atom_types)
N2_idx = find(x -> x == "N2",atom_types)
N3_idx = find(x -> x == "N3",atom_types)

V_ILJ_C_N1 = 0.0; V_ILJ_C_N2 = 0.0; V_ILJ_C_N3 = 0.0;
V_ind_C_N  = 0.0

for i in C_idx
    R1 = norm(atom_coords[i,:] - atom_coords[N1_idx[1],:])
    R2 = norm(atom_coords[i,:] - atom_coords[N2_idx[1],:])
    R3 = norm(atom_coords[i,:] - atom_coords[N3_idx[1],:])
    V_ILJ_C_N1 += ILJ_C_N1(R1)
    V_ILJ_C_N2 += ILJ_C_N2(R2)
    V_ILJ_C_N3 += ILJ_C_N1(R3)
    V_ind_C_N  += ind_C_N(R1, R2, R3)
end

# Sum everything in meV
V_ILJ_tot = (V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3)
V_int_tot = V_ILJ_tot + V_ind_C_N

# Transform in kcal/mol
V_ILJ_tot_kcal = V_ILJ_tot*0.02306
V_ind_C_N_kcal = V_ind_C_N*0.02306
V_int_tot_kcal = V_int_tot*0.02306

#println("Total Improved Lennard-Jones interaction energy: ",round(V_ILJ_tot_kcal,2)," [kcal/mol]")
@printf("\n%50s %8.2f [meV]\n","Total Improved Lennard-Jones interaction energy:",V_ILJ_tot)
@printf("%50s %8.2f [kcal/mol]\n\n","Total Improved Lennard-Jones interaction energy:",V_ILJ_tot_kcal)

#println("Induction term interaction energy: ",round(V_ind_C_N_kcal,2)," [kcal/mol]")
@printf("%50s %8.2f [meV]\n","Induction term interaction energy:",V_ind_C_N)
@printf("%50s %8.2f [kcal/mol]\n\n","Induction term interaction energy:",V_ind_C_N_kcal)

#println("Total interaction energy: ",round(V_int_tot,2)," [kcal/mol]")
@printf("%50s %8.2f [meV]\n","Total interaction energy:",V_int_tot)
@printf("%50s %8.2f [kcal/mol]\n\n","Total interaction energy:",V_int_tot_kcal)


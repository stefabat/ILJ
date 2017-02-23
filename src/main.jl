#script to execute int energy
include("ILJ.jl")
include("xyz_parser.jl")

# Define parameters
ϵ_C_N1 = 5.205; ϵ_C_N2 = 3.536;
r0_C_N1 = 3.994; r0_C_N2 = 3.818;
β_C_N1 = 8.0; β_C_N2 = 8.0;
m_C_N1 = 6.0; m_C_N2 = 6.0;
q1 = -0.56; q2 = +0.12; q3 = -0.56;
α_c = 1.2

# Define ILJ potentials
ILJ_C_N1 = ILJ_kernel(ϵ_C_N1, r0_C_N1, β_C_N1, m_C_N1)
ILJ_C_N2 = ILJ_kernel(ϵ_C_N2, r0_C_N2, β_C_N2, m_C_N2)

# Define induction potential
ind_C_N = N3_induction_kernel(r0_C_N2, β_C_N2,q1, q2, q3, α_c)

# Read geometry
atom_types,atom_coords = xyz_parser("data/n3_t555.xyz")

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

V_ILJ_tot = V_ILJ_C_N1 + V_ILJ_C_N2 + V_ILJ_C_N3
print("Total Improved Lennard-Jones interaction energy (kcal/mol): ",V_ILJ_tot*0.02306)

print("Induction term interaction energy (kcal/mol): ",V_ind_C_N*0.02306)

print("Total interaction energy (kcal/mol): ",(V_ILJ_tot+V_ind_C_N)*0.02306)

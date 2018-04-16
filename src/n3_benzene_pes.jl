
include("ILJ.jl")
include("xyz_parser.jl")

# in meV
ϵ_C_Next = 5.205; rm_C_Next = 3.994;      # Carbon to external nitrogens
ϵ_C_Nint = 3.536; rm_C_Nint = 3.818;      # Carbon to internal nitrogen
ϵ_H_Next = 2.827; rm_H_Next = 3.644;
ϵ_H_Nint = 2.431; rm_H_Nint = 3.348;
β = 8.0
m = 6.0
q_Next = -0.56; q_Nint = +0.12; # in au units
α_C = 1.20 # in Ang^3
α_H = 0.38
    
# Read the two fragments
frag_1 = readxyz("../data/geometries/benzene.xyz")
frag_2 = readxyz("../data/geometries/n3.xyz")

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

# initialize variables to save results
r = collect(linspace(3.0,8.0,501))
Vilj = zeros(length(r))
Vind = zeros(length(r))
Vels = zeros(length(r))

# generate PES
for (j,d) in enumerate(r)
    # set distance between N3- and benzene
    frag_2.coords[:,3] = d

    # Loop over fragment 1
    for i in range(1,frag_1.n_atoms)
        ri1 = vecnorm(frag_1.coords[i,:] - frag_2.coords[1,:])
        ri2 = vecnorm(frag_1.coords[i,:] - frag_2.coords[2,:])
        ri3 = vecnorm(frag_1.coords[i,:] - frag_2.coords[3,:])

        Vels[j] += els(frag_1.charges[i], frag_2.charges[1], ri1)
        Vels[j] += els(frag_1.charges[i], frag_2.charges[2], ri2)
        Vels[j] += els(frag_1.charges[i], frag_2.charges[3], ri3)
        if frag_1.types[i] == "C"
            Vilj[j] += ILJ_C_Next(ri1)
            Vilj[j] += ILJ_C_Nint(ri2)
            Vilj[j] += ILJ_C_Next(ri3)
            Vind[j] += ind_C_N3m(ri1, ri2, ri3)
        elseif frag_1.types[i] == "H"
            Vilj[j] += ILJ_H_Next(ri1)
            Vilj[j] += ILJ_H_Nint(ri2)
            Vilj[j] += ILJ_H_Next(ri3)
            Vind[j] += ind_H_N3m(ri1, ri2, ri3)
        end
    end
end


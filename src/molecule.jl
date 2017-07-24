"""
    Molecule(n_atoms, types, coords, atomic_charges)

Composite type describing a molecule

# Arguments
* `n_atoms::Integer`: Number of atoms in the molecule
* `types::Array{AbstractString}`: Nx1 array of atom symbols
* `coords::Array{Float64}`: Nx3 matrix of coordinates of the molecule
* `charges::Array{Float64}`: Nx1 array of partial atomic charges
"""
mutable struct Molecule
    n_atoms::Integer
    types  ::Array{AbstractString}
    coords ::Array{Float64}
    charges::Array{Float64}
end

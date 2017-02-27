# Simple xyz parser

function xyz_parser(filename)
    if ! isfile(filename)
        throw(ArgumentError("The file does not exist!"))
    elseif ! endswith(filename,".xyz")
        throw(ArgumentError("The file has to be a .xyz file!"))
    end
    f = open(filename,"r")
    n_atoms = parse(readline(f))
    atoms = readdlm(f,skipstart=1)
    close(f)
    if n_atoms != size(atoms,1)
        throw(ArgumentError("The number of atoms listed does not match with the number declared!"))
    end
    atom_types  = convert(Array{String,1}, atoms[:,1])
    atom_coords = convert(Array{Float64,2}, atoms[:,2:end])

    return atom_types,atom_coords
end
include("molecule.jl")

# Function to parse an xyz file
function readxyz(filename)
    if !isfile(filename)
        throw(ArgumentError("The file does not exist!"))
    elseif !endswith(filename,".xyz")
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

    return Molecule(n_atoms, atom_types, atom_coords, zeros(n_atoms,1))
end


# Check if dependency is satisfied
#if ! isdefined(:PyCall)
#  try
#    println("Loading the PyCall module...")
#    using PyCall
#  catch
#    error("You need to install the PyCall module using 'Pkg.add("PyCall")'!")
#  end
#end

using PyCall
# Import cclib to parse gaussian log file
@pyimport cclib.parser as ccp

# Function to read atomic charges from Gaussian log files
function atomic_charges(log_file, charges = "natural")
    if ! isfile(log_file)
        throw(ArgumentError("The file does not exist!"))
    elseif ! endswith(log_file,".log")
        throw(ArgumentError("The file has to be a .log file!"))
    end
    p = ccp.Gaussian(log_file)
    data = p[:parse]()
    return data[:atomcharges][charges]
end

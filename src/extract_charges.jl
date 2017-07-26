using PyCall
@pyimport cclib.parser as ccp
files = readdir()
for file in files
  if endswith(file,".log")
    name = splitext(file)[1]
    p = ccp.Gaussian(file)
    data = p[:parse]()
    writedlm(string(name,"_npa.dat"),data[:atomcharges]["natural"])
  end
end

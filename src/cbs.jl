
# Check if dependency is satisfied
if !isdefined(:LsqFit)
  try
    println("Loading the LsqFit module...")
    using LsqFit
  catch
    error("You need to install the LsqFit module!")
  end
end

# Fit SCF energies using a two- or three-point extrapolation
# formula, with optimized exponents taken from ORCA 4
function fit_scf(X,energies,basis)

  println("\nExtrapolating using the formula:\n\n\tE(X) = E(CBS) + A*exp(−αX)\n")
  
  # Two-point extrapolation
  if length(energies) == 2
    xn = [X,X+1]
    p0 = [energies[end],100.0]
    
    if X < 2
      error("X has to be at least 2!")
    
    # Starting at DZ level
    elseif X == 2
      if basis == "cc" || basis == "aug-cc"
        a = 4.42
      elseif basis == "pc"
        a = 7.02
      elseif basis == "def2"
        a = 10.39
      elseif basis == "ano"
        a = 5.41
      else
        error("Valid basis set types are:\ncc, aug-cc, pc, def2, ano")
      end
    
    # Starting a TZ level
    elseif X == 3
      if basis == "cc" || basis == "aug-cc"
        a = 5.46
      elseif basis == "pc"
        a = 9.78
      elseif basis == "def2"
        a = 7.88
      elseif basis == "ano"
        a = 4.48
      else
        error("Valid basis set types are:\ncc, aug-cc, pc, def2, ano")
      end

    else
      error("Two-point extrapolation starting at X >= 4 is not supported")
    end

  # Three-point extrapolation
  elseif length(energies) == 3
    xn = [X,X+1,X+2]
    p0 = [energies[end],100.0,5.0]
  else
    error("CBS extrapolation is supported only with two or three points")
  end
  
  if length(xn) == 2
    #model(x, p) = p[1] + p[2]*exp(-a.*sqrt(x))
    fit = curve_fit((x,p)->p[1] + p[2]*exp.(-a.*sqrt.(x)) ,xn,energies,p0)
  elseif length(xn) == 3
    #model(x, p) = p[1] + p[2]*exp(-p[3].*sqrt(x))
    fit = curve_fit((x,p)->p[1] + p[2]*exp.(-p[3].*sqrt.(x)),xn,energies,p0)
  else
    error("Contact the developer")
  end

  #println("\nExtrapolated SCF energies")
  for i=1:length(energies)
    @printf("X = %i %14s %16.8f Eh\n",X+i-1,"Energy:",energies[i])
  end
  @printf("X = CBS %12s %16.8f Eh\n","Energy:",fit.param[1])
  
  println("\nEstimated value of A: ",round(fit.param[2],2))
  if length(energies) == 3
    println("Estimated value of α: ",round(fit.param[3],2))
  end

  return fit.param[1]

end

# Fit correlation energies using a two-point extrapolation
# formula, with optimized parameters taken from ORCA 4
function fit_corr(X,energies)

  println("\nExtrapolating using the formula:\n
          E(CBS) = (X^β E(X) - Y^β E(Y)) / (X^β - Y^β)\n")
  
  if X < 2
    error("X has to be at least 2!")
  elseif X == 2
    b = 2.4     # Exponent taken from ORCA 4
  elseif X > 2
    b = 3.0     # Exponent taken from ORCA 4
  end

  E_CBS = (X^b * energies[1] - (X+1)^b * energies[2])/(X^b - (X+1)^b)

  for i=1:length(energies)
    @printf("X = %i %14s %16.8f Eh\n",X+i-1,"Energy:",energies[i])
  end
  @printf("X = CBS %12s %16.8f Eh\n","Energy:",E_CBS)
  
  return E_CBS

end

# Fit CCSD(T) energies assuming correlation error between
# MP2 and CCSD(T) is basis set independent
function fit_ccsdt(e_ccsdt, e_mp2, e_mp2_cbs)

  println("\nExtrapolating using the formula:\n
          E(CCSD(T),CBS) = E(CCSD(T),X) - E(MP2,X) + E(MP2,CBS) \n")

  E_CBS = e_ccsdt - e_mp2 + e_mp2_cbs

  @printf("X = X %14s %16.8f Eh\n","Energy:",e_ccsdt)
  @printf("X = CBS %12s %16.8f Eh\n","Energy:",E_CBS)

  return E_CBS

end

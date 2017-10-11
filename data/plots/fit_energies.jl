using LsqFit

# Define model function
model(x,p) = (p[1] .* x)./(p[2] .+ x)

# RI-SCS-MP2/cc-pVDZ interaction energies
x_mp2 = [3.0,5.0,7.0,9.0,11.0,13.0]
y_mp2 = [-26.94,-36.06,-40.81,-44.31,-47.12,-49.16]

# approximate DLPNO-CCSD(T)/cc-pVTZ interaction energies
x_cc = [3.0,5.0,7.0]
y_cc = [-36.42,-44.34,-48.94]

# guess parameters
p0=[-65.0,3.5]

# Fit model to data points
# NOTE: I get the same when fitting using gnuplot
fit_mp2 = curve_fit(model,x_mp2,y_mp2,p0)
fit_cc  = curve_fit(model,x_cc ,y_cc ,p0)

# Compute coefficient of determination as a measure
# of quality of the fit (not so sure it makes so
# much sense when having so few data points...)
ss_tot_mp2 = sum((y_mp2 - mean(y_mp2)).^2)
ss_res_mp2 = sum((fit_mp2.resid).^2)
R_mp2 = 1 - ss_res_mp2/ss_tot_mp2

ss_tot_cc = sum((y_cc - mean(y_cc)).^2)
ss_res_cc = sum((fit_cc.resid).^2)
r_cc = 1 - ss_res_cc/ss_tot_cc


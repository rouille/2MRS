PRO FitSelectionFunc, X, par, F

;;; From Crook et al. 2007
a = par[0]
b = par[1]
beta = par[2]
S = par[3] ; in Mpc
N0 = par[4]

F = N0*(beta*X/((beta*X)^b+S^b)^(1/b))^a

END
   
  
  
FUNCTION FitNGalaxiesCum, x, y, parameters = parameters, sigmaParameters = sigmaParameters, chi2 = chi2

par = [2.,4.,0.9,100.,max(y)]
fit = curvefit(x, y, fltarr(n_elements(x))+1., par, sigmaParameters, chisq = chi2, function_name = 'fitSelectionFunc', /double, /noderivative)

parameters = par

return, fit
end

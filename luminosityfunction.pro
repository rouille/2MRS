FUNCTION LuminosityFunction, Cat, magStart, magStop, magStep, distMax, magMin, mag = mag, dPhi = dPhi

;;; The luminosity function is the number of galaxies in a given
;;; volume V, of a given luminosity L. The luminosity is measured in a
;;; given waveband and thus there is a different luminosity function
;;; for every measured waveband. In a volume limited sample, the
;;; volume V is just the volume of the sample. In a flux limited
;;; sample one has to compensate for the fact that one can see
;;; brighter galaxies farther away. This can be done by determining
;;; the maximum volume that a galaxy could have been seen given its
;;; luminosity and then weighting each galaxy accordingly (cf.
;;; VolumeLimitedSample.pro with magMin and distMax).


;;; Number of step in magnitude
Nstep = abs(magStop-magStart)/magStep

;;; Distribution of the magnitude
newhist, Cat.AbsoluteMagnitude, magStart, magStop, Nstep, /noplot, xHist = mag, yHist = NGalaxies
Rmax = interpol(distMax,magMin,mag) ;;; maximum distance for AbsoluteMagnitude in Mpc
Vmax = (4./3.)*!pi*Rmax^3 ;;; associated volume in cubic Mpc

Phi = NGalaxies/Vmax/magStep
dPhi = sqrt(NGalaxies)/Vmax/magStep

return, Phi

END


FUNCTION SchechterFunc, X, alpha = alpha, MStar = MStar, PhiStar = PhiStar

;;; From Crook et al. 2007. It is a Schechter function
h = 0.71 ;;; reduced Hubble constant
alpha = -1.02
MStar = -24.1
PhiStar = 0.0108*h*h*h ;;; per cubic Mpc

F = 0.4*alog(10d)*PhiStar*10d^(0.4*(alpha+1d)*(MStar-X))*exp(-10d^(0.4*(MStar-X)))

return, F

END


FUNCTION DistanceToRedshift, distance

;;; Hubble Law in the linear regime (redshifth small regarding to 1)
;;; distance express in Mpc

Ho = 71d ;;; in km/s/Mpc according to WMAP seven years
c = 299792458d ;;; speed of light in m/s

velocity = Ho*distance ;;; recessional velocity expressed in km/s.

redshift = velocity/(c/1000d)

return, redshift
END



FUNCTION RedshiftToDistance, redshift

;;; Hubble Law in the linear regime (redshifth small regarding to 1)
;;; distance express in Mpc

Ho = 71d ;;; in km/s/Mpc according to WMAP seven years
c = 299792458d ;;; speed of light in m/s

velocity = redshift*(c/1000d) ;;; recessional velocity expressed in km/s.

distance = velocity/Ho

return, distance
END



FUNCTION ApparentMagnitudeToFlux, ApparentMagnitude

;;; Magnitude to flux density converter

;;; 2MASS isophotal flux-for-0-magnitude from Cohen et al. (2003)
F0 = 666.8 ; in Jy
Flux = F0*10^(-ApparentMagnitude/2.5) ; in Jy

return, flux
END



FUNCTION AbsoluteMagnitudeToApparentMagnitude, AbsoluteMagnitude, distance

ApparentMagnitude = AbsoluteMagnitude+25.+5.*alog10(distance)

return, ApparentMagnitude
END



FUNCTION ApparentMagnitudeToAbsoluteMagnitude, ApparentMagnitude, distance

AbsoluteMagnitude = ApparentMagnitude-25.-5.*alog10(distance)

return, AbsoluteMagnitude
END



FUNCTION DistanceModulusToDistance, magnitude

distance = 10^(0.2*magnitude+1) ;;; in pc
distance = distance/1d6

return, distance
END



FUNCTION AbsoluteMagnitudeToLuminosity, AbsoluteMagnitude

AbsoluteMagnitudeSun = 3.28 ;;; absolute magnitude of the sun in the K band
Luminosity = 10d^((AbsoluteMagnitudeSun-AbsoluteMagnitude)/2.5) ;;; luminosity in suns

return, Luminosity
END

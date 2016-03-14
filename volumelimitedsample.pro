PRO VolumeLimitedSampleData, Cat, distStart, distStop, distStep, distMax = distMax, magMin = magMin, size = size

;;; Return three tables:
;;; Galaxies with absolute magnitude greater than magMin will define a
;;; volume limited sample up to distance equal to distMax. The table size
;;; gives the the number of galaxies in the volume limited sample
;;; defined by [magMin,distMax]


;;; Number of step in redshift
Nstep = abs(distStop-distStart)/distStep


;;; Output tables
magMin = fltarr(Nstep)
distMax = fltarr(Nstep)
size = fltarr(Nstep)

;;; Filling in the tables
for i = 0l, Nstep-1 do begin
    wDist = where(Cat.distance ge distStart+i*distStep and Cat.distance lt distStart+(i+1)*distStep)
    if wDist[0] ne -1 then begin
        wMagMax = where(Cat[wDist].AbsoluteMagnitude eq max(Cat[wDist].AbsoluteMagnitude))
        magMin[i] = Cat[wDist[wMagMax[0]]].AbsoluteMagnitude ;;; several entries with same magMin is possible
        distMax[i] = max(Cat[wDist[wMagMax]].distance)
    endif
endfor

for i = 0l, Nstep-1 do size[i] = n_elements(where(Cat.distance lt distMax[i] and Cat.AbsoluteMagnitude lt magMin[i]))


END





PRO VolumeLimitedSampleAnalytic, Cat, distStart, distStop, distStep, distMax = distMax, magMin = magMin, size = size, m0 = m0

;;; Return three tables:
;;; Galaxies with absolute magnitude greater than magMin will define a
;;; volume limited sample up to distance equal to distMax. The table size
;;; gives the the number of galaxies in the volume limited sample
;;; defined by [magMin,distMax]

;;; 2RMS flux threshold
if not keyword_set(m0) then m0 = 11.25


;;; Number of step in redshift
Nstep = abs(distStop-distStart)/distStep


;;; Output tables
magMin = fltarr(Nstep)
distMax = fltarr(Nstep)
size = fltarr(Nstep)

for i = 0l, NStep-1 do distMax[i] = 0.1+distStart+i*distStep
for i = 0l, NStep-1 do MagMin[i] = ApparentMagnitudeToAbsoluteMagnitude(m0,distMax[i])
for i = 0l, Nstep-1 do size[i] = n_elements(where(Cat.distance lt distMax[i] and Cat.AbsoluteMagnitude lt magMin[i]))

END


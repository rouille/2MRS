FUNCTION read_2MRS, FileName

print, 'Reading '+FileName
readcol, FileName, PGC, ApparentMagnitude, velocity, LonGal, LatGal, Luminosity, Boost, LuminosityBoosted, DistanceCorrected, skipline = 1, delimiter = '|'
keep = where(ApparentMagnitude ne 0 and velocity ne 0 and LonGal ne 0 and LatGal ne 0 and Luminosity ge 1d7)


;;; Structure
CatTemplate = {ra:0d, dec:0d, l:0d, b:0d, distance:0d, Flux:0d, ApparentMagnitude:0d, AbsoluteMagnitude:0d, Luminosity:0d}
Cat = replicate(CatTemplate,n_elements(keep))

;;; Speed of light
c = 299792458d ;;; speed of light in m/s

;;; Fill out the structure
for i = 0l, n_elements(keep)-1 do begin
    Cat[i].l = LonGal[keep[i]]
    Cat[i].b = LatGal[keep[i]]
    euler, LonGal[keep[i]], LatGal[keep[i]], ra, dec, 2
    Cat[i].ra = ra
    Cat[i].dec = dec
    if DistanceCorrected[keep[i]] ne 0 then begin
        distance = DistanceCorrected[keep[i]]
        Cat[i].Distance = distance
    endif else begin
        distance = RedshiftToDistance(velocity[keep[i]]*1000d/c)
        Cat[i].Distance = distance
    endelse
    Cat[i].ApparentMagnitude = ApparentMagnitude[keep[i]]
    Cat[i].Flux = ApparentMagnitudeToFlux(ApparentMagnitude[keep[i]])
    Cat[i].AbsoluteMagnitude = ApparentMagnitudeToAbsoluteMagnitude(ApparentMagnitude[keep[i]],distance)
    Cat[i].Luminosity = Luminosity[keep[i]] ;;; in suns
endfor


return, Cat
END

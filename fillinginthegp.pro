FUNCTION FillingInTheGP, Cat, distMaxObs, MagMinObs

;;; According to Crook et al. 2007



;;; Structure
CatTemplate = {ra:0d, dec:0d, l:0d, b:0d, distance:0d, Flux:0d, ApparentMagnitude:0d, AbsoluteMagnitude:0d, Luminosity:0d}



;;; Reduced Hubble constant: Ho = 100 * h km/s/Mpc
h = 0.71 ; according to WMAP seven years
;;; Speed of light in m/s
SpeedOfLight = 299792458*1d
;;; Scatter in redshift/distance for the generated galaxies
RedshiftScatter = 20./(1d*SpeedOfLight/1000)
DistScatter = RedshiftToDistance(RedshiftScatter)



;;; We divide the catalog into bins spanning 10° in longitude and 10/h in distance
nLon = 36
LonBin = 360/36
DistBoxSize = 10.
nDistance = ceil(max(Cat.distance)/(DistBoxSize/h))


FirstGalaxy = 'true'
print, '###############################'
for i = 0l, nLon-1 do begin
    LonMin = i*10
    LonMax = (i+1)*10
    if LonMax le 30 or LonMin ge 330 then begin
        LatMin = 10
        LatMax = 20
        factor = 1
    endif else begin
        LatMin = 5
        LatMax = 15
        factor = 0.5
    endelse
    print, strtrim(LonMin,2)+' < l < '+strtrim(LonMax,2)+' and '+strtrim(LatMin,2)+' < |b| < '+strtrim(LatMax,2)
    for j = 0l, nDistance-1 do begin
        DistMin = DistBoxSize*j/h
        DistMax = DistBoxSize*(j+1)/h
        print, '+++ '+ strtrim(DistMin,2)+' <= d(Mpc) < '+strtrim(DistMax,2)+' +++'
        wAdjacent = where(Cat.l ge LonMin and Cat.l lt LonMax and abs(Cat.b) gt LatMin and abs(Cat.b) lt LatMax and Cat.distance ge DistMin and Cat.distance lt DistMax)
        wIn = where(Cat.l ge LonMin and Cat.l lt LonMax and abs(Cat.b) le LatMin and Cat.distance ge DistMin and Cat.distance lt DistMax)
        if wAdjacent[0] ne -1 then nAdjacent = n_elements(wAdjacent) else nAdjacent = 0
        if wIn[0] ne -1 then nIn = n_elements(wIn) else nIn = 0
        print, 'Number of Galaxies in the GP: '+strtrim(nIn,2)
        print, 'Number of Galaxies in the adjacent bins: '+strtrim(nAdjacent,2)+' (factor = '+strtrim(factor,2)+')'
        nDrawn = 0 ;;; Number of galaxies drawn in the bin
        if nAdjacent ge 1 then begin
            nDrawn = nAdjacent*factor+randomu(seed,1,/normal)
            nDrawn = long(abs(nDrawn[0])) ;;; it seems that in Crook et al. 2007, nIn is not substracted
            if nDrawn ge 1 then begin
                CatGeneratedTmp = replicate(CatTemplate,nDrawn)
                CatGeneratedTmp.b = 2*LatMin*randomu(seed, nDrawn, /uniform)-LatMin ;;; Latitude (random in the plane)
                for k = 0l, nDrawn-1 do begin
                    RandomIndex = RandPerm(nAdjacent)
                    CatGeneratedTmp[k].l = Cat[wAdjacent[RandomIndex[0]]].l ;;; Longitude (pick one in the adjacent bins)
                    distanceTmp = Cat[wAdjacent[RandomIndex[0]]].distance ;;; distance (pick one in the adjacent bin)
                    CatGeneratedTmp[k].distance = distanceTmp+DistScatter*randomu(seed,1,/normal) ;;; additional scatter for the distance
                endfor
                if FirstGalaxy eq 'true' then begin 
                    CatGenerated = CatGeneratedTmp
                    FirstGalaxy = 'false'
                endif else begin
                    CatGenerated = [CatGenerated,CatGeneratedTmp]
                endelse
            endif else begin
                nDrawn = 0
            endelse
        endif
        print, 'Number of Galaxies drawn: '+strtrim(nDrawn,2)
    endfor
    print, '###############################'
endfor


;;; Equatorial coordinates
euler, CatGenerated.l, CatGenerated.b, ra, dec, 2
CatGenerated.ra = ra ;;; Right Ascension
CatGenerated.dec = dec ;;; Declination



;;; Absolute Magnitude
magMinSource = interpol(magMinObs,distMaxObs,CatGenerated.distance) ;;; minimum magnitude for the source
for i = 0l, n_elements(magMinSource)-1 do begin
    wKeep = where(Cat.AbsoluteMagnitude le magMinSource[i] and Cat.distance le CatGenerated[i].distance)
    if wKeep[0] ne -1 then begin
        RandomIndex = RandPerm(n_elements(wKeep))
        CatGenerated[i].AbsoluteMagnitude = Cat[wKeep[RandomIndex[0]]].AbsoluteMagnitude ;;; Absolute Magnitude
        CatGenerated[i].Luminosity = Cat[wKeep[RandomIndex[0]]].Luminosity ;;; Luminosity in Suns
    endif
endfor



;;; Apparent Magnitude
CatGenerated.ApparentMagnitude = AbsoluteMagnitudeToApparentMagnitude(CatGenerated.AbsoluteMagnitude,CatGenerated.distance)


;;; Flux
CatGenerated.Flux = ApparentMagnitudeToFlux(CatGenerated.ApparentMagnitude)

return, CatGenerated

END





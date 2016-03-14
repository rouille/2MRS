FUNCTION Complete, Cat, distRef, magRef, distComplete

;;; Assuming that Cat is a volume-limited-sample of sources up to
;;; distRef, this code generates a catalog with sources homogeneously
;;; distributed on the sky in the range [distRef,distComplete] (in Mpc). 
;;; The number of sources is set according to the density of sources in
;;; the original sample. The magnitude of the sources generated is
;;; drawn randomly among the distribution of magnitude of the original
;;; catalog.



;;; Structure
CatTemplate = {ra:0d, dec:0d, l:0d, b:0d, distance:0d, Flux:0d, ApparentMagnitude:0d, AbsoluteMagnitude:0d, Luminosity:0d}



;;; density of the original sample
VolumeRef = (4./3.)*!pi*(1d*distRef)^3
density = n_elements(Cat)/VolumeRef



;;; Number of sources to draw between [distRef,distForComplete]
VolumeComplete = (4./3.)*!pi*(1d*distComplete)^3
nSourcesInTotalVolume = round(density*VolumeComplete)
nSourcesToDraw = nSourcesInTotalVolume-n_elements(Cat)
print, strtrim(nSourcesToDraw)+'  will be drawn'


;;; Current number of sources drawn
nSourcesCurrent = 0l
CatGenerated = replicate(CatTemplate,nSourcesToDraw)


while nSourcesCurrent lt nSourcesToDraw do begin
    uvTmp = distComplete*(2*randomu(seed,3)-1)
    rTmp = sqrt(uvTmp[0]*uvTmp[0]+uvTmp[1]*uvTmp[1]+uvTmp[2]*uvTmp[2])
    if rTmp le distComplete and rTmp ge distRef then begin
        x = uvTmp[0]
        y = uvTmp[1]
        z = uvTmp[2]
        theta = acos(z/rTmp)
        if x ne 0 then begin
            phi = atan(y,x)
        endif else begin
            phi = !pi/2d
        endelse
        if phi lt 0 then phi = phi +2d*!pi

        CatGenerated[nSourcesCurrent].distance = rTmp
        CatGenerated[nSourcesCurrent].l = phi*1d/!dtor
        CatGenerated[nSourcesCurrent].b = 90.-theta*1./!dtor
        nSourcesCurrent = nSourcesCurrent+1l
        if 50000*(nSourcesCurrent/50000) eq nSourcesCurrent then print, '# '+strtrim(nSourcesCurrent,2)+' sources drawn'
    endif
endwhile



;;; Equatorial coordinates
euler, CatGenerated.l, CatGenerated.b, CatGenerated.ra, CatGenerated.dec, 2



;;; Magnitude and flux
magMax = floor(min(Cat.AbsoluteMagnitude))

mag = magRef+(magMax-magRef)*findgen(100)/99d

newhist, Cat.AbsoluteMagnitude, magMax, magRef, 20, xhist = x, yhist = y, /noplot
w = where(y ne 0)
phi = interpol(y[w],x[w],mag)

AbsoluteMagnitudeTmp =  distri(mag,phi,nSourcesToDraw*2l)
w = where(AbsoluteMagnitudeTmp lt magRef)
AbsoluteMagnitude = fltarr(nSourcesToDraw)
for i = 0l, nSourcesToDraw-1 do AbsoluteMagnitude[i] = AbsoluteMagnitudeTmp[w[i]]

CatGenerated.AbsoluteMagnitude = AbsoluteMagnitude
CatGenerated.Luminosity = AbsoluteMagnitudeToLuminosity(CatGenerated.AbsoluteMagnitude)
CatGenerated.ApparentMagnitude = AbsoluteMagnitudeToApparentMagnitude(CatGenerated.AbsoluteMagnitude,CatGenerated.distance)
CatGenerated.Flux = ApparentMagnitudeToFlux(CatGenerated.ApparentMagnitude)


return, CatGenerated
END





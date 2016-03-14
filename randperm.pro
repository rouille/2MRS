FUNCTION RandPerm, numberOfElements, seed = seed

x = lindgen(numberOfElements)

index = x[Sort(Randomu(seed, numberOfElements))]

return, index
END

# Mass Balance Calculations 
# Shawn.Salley@usda.gov, 20230918

# Strain is the amount or sense of deformation in soils using an assumed isovolume frame 
# of reference, such as an insoluble host mineral of immobile element (i, usually Titanium
# or Zirconium). Positive strains infer dilation and negative strains represent collapse.

# Strain Calculation: using the Strain (ε) formula # 2, Brimhall et al., 1991
# variables = immobile element (ie), bulk density (p), parent material (pm), weathered horizon (w)
# formula: strain = [(pm.p)*(pm.ie)] / [(w.p)*(w.ie)]

# test data

library(aqp)
#
hz <- read.csv("D:/FY2023/Projects/Fraiser Soil/bergstrom2019_hz.csv")
depths(hz) <- Site ~ top + bottom

# function with coded columns (where deepest horizon is considered 'parent material')

f <- function(x) { ((tail(x$BD, 1) * tail(x$Ti, 1)) / (x$BD * x$Ti)) - 1 }

profileApply(hz, FUN=f)  # this matches the εTi calculation in 'bergstrom2019'


# AQP function (I need help here), it should look something like this... 

horizon_strain <- function(spc, bd, im){
  #
  f <- function(x) { ((tail(x$bd, 1) * tail(x$im, 1)) / (x$bd * x$im)) - 1 }
  #
  i.strain <- profileApply(spc, FUN=f)
  return(i.strain)
  #
}  

horizon_strain(hz, BD, Ti)

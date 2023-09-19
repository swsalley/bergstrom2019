
# Development of Soil Mass Balance Calculations ####

# Functions: Strain, MassTransfer, MassFlux, weathered.

# Shawn.Salley@usda.gov, 20230918

# Initial Testing of Mass Balance Formulas and Function Development using Bergstrom2019 ####

library(aqp)
bergstrom2019 <- read.csv("D:/FY2023/Projects/Fraiser Soil/bergstrom2019_hz.csv")
depths(bergstrom2019) <- Site ~ top + bottom



# Strain ####

# Strain (ε) using formula # 2 from Brimhall et al., 1991
# Strain is the amount or sense of deformation in soils using an assumed isovolume 
# frame of reference, such as an insoluble host mineral of immobile element (usually Titanium
# or Zirconium). Positive strains infer dilation and negative strains represent collapse.

# Variables = immobile element (ie), bulk density (p), parent material (pm), weathered horizon (w)
# formula: strain = [(pm.p)*(pm.ie)] / [(w.p)*(w.ie)]
Strain <- function(x) { ((tail(x$BD, 1) * tail(x$Ti, 1)) / (x$BD * x$Ti)) - 1 }

# Apply strain for each profile
bergstrom2019$strain <- profileApply(bergstrom2019, FUN=Strain)  

# compare hz$strain with εTi in 'bergstrom2019'
plot ( x = bergstrom2019$strain, y = bergstrom2019$X.Ti)
text ( x = bergstrom2019$strain, y = bergstrom2019$X.Ti, labels = bergstrom2019$hzID)
#  Note: Soil horizons IL5 (137:140) are off the 1:1 line (landscape position not used in 2019 paper)

# end Strain #


# element mobility (Ti or Zr) ####
# Ti remains enriched in the fine fraction and Zr enriched in the coarse silt fraction 
# of the particle size distribution (Stiles et al., 2003; Taboada et al., 2006).
# test of mobility, regress mean strain for Ti against Zr

# regressing the relationship between their mass transfer function values against percent clay and sand

# end test of mobility




# Mass Transfer ####

# Coefficient of mass transfer (T j,w) using formula 4 from Brimhall et al., 1991.
# Mass Transfer is the horizon's element mobility in the soil per mass fraction of the parent material.
# Variables = mobile element (me), bulk density (p), parent material (pm), weathered horizon (w), horizon strain (strain)
# Formula : MassTransfer = (bd.w * me.w) / (bd.pm * me.w) * (strain + 1) - 1

MassTransfer<- function(x) { ( x$BD * x$Ca ) / (tail(x$BD, 1) * tail(x$Ca, 1)) * (x$strain + 1) - 1 }

# Apply MassTransfer for each profile

bergstrom2019$T_Ca <- profileApply(bergstrom2019, FUN=MassTransfer)  

# compare hz$T_Ca with X.ca from 'bergstrom2019'
plot ( x =bergstrom2019$T_Ca,  y = bergstrom2019$X.ca)
text ( x =bergstrom2019$T_Ca,  y = bergstrom2019$X.ca, labels = bergstrom2019$hzID)

# errors ... will have to unpack later, maybe change to concentration?

# end mass transfer #




# Mass Flux ####
#
# Mass flux (mass j, flux) using formula 3 from Bergstrom et al., 2019 (from Egli etal 2000)
# Mass flux is the horizon weathering by thicknes of a mobile element through the soil profile

# Variables= profile weathering depth (wd), horizon thickness (z), strain (strain), 
#            horizon mass transfer (MassTransfer)

# Formula : MassFlux = ( bd.pm * z ) * ( 1 / (strain.w + 1)) * (mobile element.w * MassTransfer.w)




# Still To Do when get full dataset from Robert:



# 1) calculate weatherable profile depth, something like:

h <- function(x) { (tail(x$top, 1))}
hz$weathered.depth <-  profileApply(spc, FUN = h) # weathered profile depth (depth to C)
# spc.top <- trunc(spc, 0, weathered.depth) # may not need to subset

# 1 b)  thickens criteria , something like this:

hz$thick <- hz$bottom - hz$top
h <- function(x) { ( x$thick / x$weathered.depth ) } 
profileApply(hz, FUN = h)
#


spc <- hz
weathered.depth <- function(x) { (tail(x$top, 1))}
spc$weathered.depth <-  profileApply(spc, FUN = weathered.depth)  
spc$thick <- spc$bottom - spc$top
hz.weathered.pct <- function(x) { ( x$thick / x$weathered.depth ) } 
spc$hz.weathered.pct <- profileApply(spc, FUN = hz.weathered.pct) 
  

hz$thick
hz$weathered.depth

# 2) Horizon Mass Flux = horizon percent of weathered.depth (z%)
MassFlux <- function(x) { (bd * z%) * (1/(strain + 1) * (x$T_Ca ) }



# 4) Apply MassFlux for each horizon
hz$mf_Ca <-  profileApply(hz, FUN = MassFlux)


# end MassTransfer #

# end #






#graveyard
# Conversion to depth weighted mass transfer, multiply each value by percent of profile
h <- function(x) { (tail(x$top, 1))}
topsoil <-  profileApply(hz, FUN = h)  
hz.top <- trunc(hz, 0, topsoil)
# might be a better way to calculate percent of weatherable profile
hz.top$thick <- hz.top$bottom - hz.top$top
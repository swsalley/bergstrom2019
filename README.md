# Soil Mass Balance Data from Bergstrom et al. 2019.

A new SPC-based example dataset called ["bergstrom2019"](https://github.com/ncss-tech/aqp/issues/298) of soil and landscape constituent mass balance data for Fraser Experimental Forest. Dataset is to develop and test new methods for working with strain calculations and vizualizations. Data are an expanded datset that was used in [Bergstrom et al, 2019](https://www.sciencedirect.com/science/article/pii/S0016706117314738) to determine whether landscape position along catenas imparts a control on the distribution and sourcing of soil cations in the FEF and evaluate the contribution of atmospherically-derived Ca to the soil cation pool in FEF soils. A select dataset was previously published with the 2019 paper with Mendeley at [doi:10.17632/yj4jh8sv3y.1](https://doi.org/10.17632/yj4jh8sv3y.1).

## Overview

The Fraser Experimental Forest (FEF), Grand County, Colorado, USA, is located in the central Rocky Mountains (Fig. 1). Research in the fields of hydrology and forest dynamics has taken place in the FEF since 1937. 

![Fraser Site Map](https://github.com/swsalley/bergstrom2019/blob/main/map.jpg?raw=true)

Fig. 1,  Soil and catchment map of the [Fraser Experimental Forest](https://www.fs.usda.gov/main/fraser/home) near Fraser, Colorado. Data were used in calculating soil strain, mass transfer function, and atmospheric deposition contributions to soil catenas in the Fraser Experimental Forest. Elemental and physical data are presented by horizon, along with horizon thicknesses and nomenclature. Data are organized by site. Each site (BL1, BL2, etc.) is a location along a soil catena. Catena name corresponds to the fist two letters of the site name (BL, BU, EL, EU, FL, FU, IL, and IU). Landscape position along catena corresponds to the number after the site name (1=summit, 2=shoulder, 3=backslope, 4=footslope). 

## Geochemical (constituent) Mass Balance

The geochemical (constituent) mass balance approach is used to estimate weathering by calculating volume changes through a soil profile and parent material composition. The analytical functions are based on the principle of conservation of mass and include a term quantifying mass flux into/out of the soil and between horizons. More information on the approach and calculations are found in the following manuscripts (Brimhall and Dietrich, 1987; Chadwick et al., 1990; Brimhall et al., 1992). The approach uses Strain (ε), a measure of soil volume change incurred during pedogenesis. Strain index is unitless and calculated as the sum of the depth-weighted contibutions from each weathered soil horizon in its respective pedon. The mass transfer coefficient (τ) is used to evaluate element mobility within the soil, and mass flux is used to evaluated mbility within the landscape. 

## References

- Bergstrom, R. M., Borch, T., Martin, P. H., Melzer, S., Rhoades, C. C., Salley, S. W., & Kelly, E. F. (2019). The generation and redistribution of soil cations in high elevation catenas in the Fraser Experimental Forest, Colorado, US. Geoderma, 333, 135-144. [Link](https://doi.org/10.1016/j.geoderma.2018.07.024).
- Bergstrom, R. M., Borch, T., Martin, P. H., Melzer, S., Rhoades, C. C., Salley, S. W., & Kelly, E. F. (2019). Soil Data associated with manuscript "The generation and redistribution of soil cations in high elevation catenary sequences in the Fraser Experimental Forest, Colorado, U.S." In Geoderma. Mendeley. [Link](https://doi.org/10.17632/yj4jh8sv3y.1).
- Brimhall, G. H., & Dietrich, W. E. (1987). Constitutive mass balance relations between chemical composition, volume, density, porosity, and strain in metasomatic hydrochemical systems: results on weathering and pedogenesis. Geochimica et Cosmochimica Acta, 51(3), 567-587. [Link]( https://doi.org/10.1016/0016-7037(87)90070-6).
- Brimhall, G. H., Chadwick, O. A., Lewis, C. J., Compston, W., Williams, I. S., Danti, K. J., ... & Bratt, J. (1992). Deformational mass transport and invasive processes in soil evolution. Science, 255(5045), 695-702. [Link](https://doi.org/10.1126/science.255.5045.695).
- Chadwick, O. A., Brimhall, G. H., & Hendricks, D. M. (1990). From a black to a gray box—a mass balance interpretation of pedogenesis. Geomorphology, 3(3-4), 369-390. [Link](https://doi.org/10.1016/0169-555X(90)90012-F).

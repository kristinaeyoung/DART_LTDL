# DART-updated

## Description
This repository contains a version of the Disturbance Automated Reference Toolset (DART) that has been updated to use the R packages sf (cite) and terra (cite) rather than the depreciated sp package and the raster package. DART was developed by the US Geological Survey to assess ecological recovery from land disturbances such as the construction of oil and gas well pads. DART does this by matching pixels within the disturbance areas with those outside of it using edaphic and topographic variables and Gower distance (Nauman et al. 2017). The version of DART provided here matches each pixel in the disturbance polygon with candidate pixels in an iterative fashion (Fick et al. 2022). Also included is code to generate a counterfactual for each disturbed pixel using synthetic control, which provides an estimate of the response variable(s) in absence of the disturbance (Abadie et al. 2010, Fick et al. 2021). 

## Spatial Data Notes
All data provided in this repository are projected in Albers Equal Area Conic. Polygon data shapefiles provided by the user should be projected in the same CRS. The shapefile should include the following attributes: 'polyID' (a unique identifier for that polygon, preferably an integer), 'YOI' (Year of intervention in the date format YYYY), and DOI (Date of intervention, in the date format YYYY-MM-DD).

## References and further reading
Abadie, A., A. Diamond, and J. Hainmueller. 2010. Synthetic control methods for comparative case studies: Estimating the effect of California’s tobacco control program. Journal of the American statistical Association 105:493–505.

Nauman, T.W., M.C. Duniway, M.L. Villarreal, and T.B. Poitras. 2017. Disturbance automated reference toolset (DART): Assessing patterns in ecological recovery from energy development on the Colorado Plateau. Science of the Total Environment, 584, 476-488.

Nauman, T.W., and M.C. Duniway. 2016. The automated reference toolset: a soil-geomorphic ecological potential matching algorithm. Soil Science Society of America Journal 80:1317–1328.

Fick, S.E., T.W. Nauman, C.C. Brungard, and M.C. Duniway, 2021. Evaluating natural experiments in ecology: using synthetic controls in assessments of remotely sensed land treatments. Ecolological Applications 31:3 https://doi.org/10.1002/eap.v31.310.1002/eap.2264.

Fick, S.E., T.W. Nauman, C.C. Brungard, and M.C. Duniway. 2022. What determines the effectiveness of 
Pinyon-Juniper clearing treatments? Evidence from the remote sensing archive and counter-factual scenarios. Forest Ecology and Management, 505:119879.

USGS Info Page: https://www.usgs.gov/publications/disturbance-automated-reference-toolset-dart-assessing-patterns-ecological-recovery

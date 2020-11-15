# National-Anuran-Vulnerability-Assessment-Project
Based upon code from the [National Fishes Vulnerability Assessment Project](https://github.com/silknets/National-RCS)

===============================

Sam Silknetter (SCS) Updated 04/21/2020
Traci P. DuBose Updated 11/15/2020

-------------------------------

# Purpose

This repository contains R code for analyzing species rarity and climate sensitivity for 92 anuran species native to the contiguous (Lower 48) United States. We use downloaded observations from Global Biodiversity Information Facility (GBIF) and HerpMapper to describe where these species occur. Once filtered, occurrence points are then used to determine area of occupancy (AOO) for two grain sizes: small occupied watersheds (USGS HUC-12 sub-basins and HydroBASINS level 12) and buffered occurrence points (1km radius). The PRISM climate model (AN81m)[1] is then used to determine species climate sensitivity (CS) for each species' grain-specific AOO within the contiguous US. Outside of the contiguous US, we use WorldClim data to describe the climate sensitivity of each species [2]. Climate variables examined include annual precipitation (Ppt), annual minimum temperature (Tmin) and annual maximum temperature (Tmax). At each grain size, CS and AOO metrics are combined into a relative rarity and climate sensitivity index (RCS)[3].

-------------------------------

# Scope

This National Anuran Vulnerability Assessment Project provides a relative assessment of climate vulnerability for most of the native frogs and toads across the contiguous United States. We also assess how well current conservation status aligns with this measure of vulnerability.

-------------------------------

# Intended Uses

The analyses herein are intended for use by wildlife managers and conservation practitioners to identify species that have high relative, intrinsic sensitivity to changes in climate. By using relative metrics that can be applied across different spatial scales, these assessments of geographic sensitivity allow for direct comparisons between species with variable data availability, including species that are both well- and poorly-studied. One note: as a relative metric, the ranking and index value are sensitive to the spatial extent considered and the taxa included in the analysis. 

-------------------------------

# References

[1] PRISM Datasets - http://www.prism.oregonstate.edu/documents/PRISM_datasets.pdf

[2] WorldClim Datasets - https://www.worldclim.org/data/worldclim21.html 

[3] Mims, M. C., D. H. Olson, D. S. Pilliod, and J. B. Dunham. 2018. Functional and geographic components of risk for climate sensitive vertebrates in the Pacific Northwest, USA. Biological Conservation 228:183-194.

===============================

# RCS Script Workflow

-------------------------------

Step 1: Download Occurrence data

This can be done manually or through [Sam's species filtering script](https://github.com/silknets/National-RCS/blob/master/Species%20Filtering%20Code_SCS.R) in R. Chloe Moore and Tess Alexander downloaded the Anuran occurrences in September 2019 for the main analysis. The script "LoadingOccDataframes.R" can be used to generate one large csv file of all occurrence points and to check the taxonomy against ITIS. 

-------------------------------

Script 1: "Creating HUC12 RDS.R"

Code to load the USGS HUC12 spatial data from the National Watershed Database. You can download the USGS National Watershed Database [geodatabase here.](https://www.usgs.gov/core-science-systems/ngp/national-hydrography/access-national-hydrography-products) Within this script, we also create the subset of HydroBasin level 12s used to descrbe the area outside of the Contiguous US. The HydroBASINS geodatabase can be [found here.](https://hydrosheds.org/page/hydrobasins) This script also removes HUCs outside of the contiguous US and removes some artifacts from the Great Lakes region. In this analysis, estuary HUCS are still within the spatial dataframes. Outputs three RDS files: a xxx projected map of the contiguous US,  a xxx projected map of all HUC12s within the contiguous US, and a xxx projected map of the HydroBasins within the North American continent but outside the contiguous US. 

-------------------------------

Script 2: "Area of Occupancy Code_Anurans.R"

Code to calculate area of occupancy (AOO) at two grains (small watersheds and 1km Point Buffer) for all species in the final list. Output files include  shapefiles at both grains/species, and a CSV with AOO values and ranks for both grain sizes. Contains code for Figure

-------------------------------

Script 3: "ARC_PRISM Download_SCS.R"

This code downloads PRISM data using library(prism). Data may be retrieved through the PRISM website - see https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf for additional details.  

-------------------------------

Script 4: "Climate Sensitivity_Anurans.R"

This code must be run on the super computer (need higher working memory). Code to calculate standard deviations of climate variables for each species. This requires PRISM rasters downloaded via Script 3 and WorldClim rasters as well. 

-------------------------------

Script 5: "RCS_Index_Code.R"

This code calculates Climate Sensitivity (CS) from standard deviation values generated in Script 4. Area of Occupancy (AOO) data (generated in Script 2) and CS data are merged into a single RCS output for all species and grains. Finally, the Relative Climate Sensitivity index (RCS) is calculated for each species and range metric. 

-------------------------------

Script 6: "RCS_Conserv_Status.R"

This code calculates Climate Sensitivity (CS) from standard deviation values generated in Script 4. Area of Occupancy (AOO) data (generated in Script 2) and CS data are merged into a single RCS output for all species and grains. Finally, the Relative Climate Sensitivity index (RCS) is calculated for each species and range metric. 

===============================

# Other Files:

-------------------------------
RCS_decision_notes.Rmd

This is a document that describes some early decisions when starting this project. It discusses why some taxa are excluded (those with < 20% of their AOO within the contiguous US), why we mostly consider the entire extent, native range extent for the anurans, and a brief comparison of the WorldClim and PRISM datasets. We also explore how the occurrence data shifts through time and how individuals can have marked impacts on certain taxa. 

-------------------------------
make_occ_maps.Rmd

This creates a document that has the IUCN range maps and occurrence points for each anuran species considered within this analysis. 

-------------------------------
anuran_rcs_qaqc.Rmd 

This file documents the quality assurrence and quality control measures used to ensure the precision of the RCS calculated using the workflow above. It will be Appendix 3 in the coming manuscript. It should only be run after completing the workflow above as it relys on outputs from Script 1, Script 2, and Script 4

-------------------------------

TODO: add input and output files to the above scripts.
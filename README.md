# National-Anuran-Vulnerability-Assessment-Project
Based upon code from the [National Fishes Vulnerability Assessment Project](https://github.com/silknets/National-RCS)

===============================

Sam Silknetter (SCS) Updated 04/21/2020
Traci P. DuBose Updated 04/22/2021

-------------------------------

# Purpose

This repository contains R code for analyzing species rarity and climate sensitivity for 90 anuran species native to the contiguous (Lower 48) United States. We use downloaded observations from Global Biodiversity Information Facility (GBIF) and HerpMapper to describe where these species occur. Once filtered, occurrence points are then used to determine area of occupancy (AOO) for two grain sizes: small occupied watersheds (USGS HUC-12 sub-basins and HydroBASINS level 12) and buffered occurrence points (1km radius). The PRISM climate model (AN81m)[1] is then used to determine species climate sensitivity (CS) for each species' grain-specific AOO within the contiguous US. Outside of the contiguous US, we use WorldClim data to describe the climate sensitivity of each species [2]. Climate variables examined include annual precipitation (Ppt), annual minimum temperature (Tmin) and annual maximum temperature (Tmax). At each grain size, CS and AOO metrics are combined into a relative rarity and climate sensitivity index (RCS)[3].

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

--------------------------------------------------------------------------

Build the spatial datasets at the designated grain size -- "0_Creating HUC12 RDS.R"

Code that creates the spatial rds files used to quantify area of occupancy. You can download the USGS National Watershed Database [geodatabase here.](https://www.usgs.gov/core-science-systems/ngp/national-hydrography/access-national-hydrography-products) Within this script, we also create the subset of HydroBasin level 12s used to descrbe the area outside of the Contiguous US. The HydroBASINS geodatabase can be [found here.](https://hydrosheds.org/page/hydrobasins) In order to still incorporate USGS HUC-12 watershed subbasins, we stiched together HydroBasins's level 12 watershed polygons from outside the contiguous US with USGS's HUC-12 polygons. Each polygon has a unique ID that relates it to the larger watershed area. This script 1) isolates USGS watersheds within the contiguous US, 2) removes overlapping polygons between the HydroBasins and HUC-12 datasets, and 3) creates HydroBasin clipped polygons that bridge the two datasets described above. It also outputs three RDS files: a projected map of the contiguous US, a projected map of all HUC12s within the contiguous US, and a projected map of HydroBasins within the North American continent but outside the contiguous US and our created clipped HydroBasin polygons.

-------------------------------

Download Occurrence data -- "0_GBIF_occurrence_download.R"

We downloaded GBIF occurrences using the rgbif R package to query the GBIF API. Chloe Moore and Tess Alexander downloaded the HerpMapper occurrences in September 2019. These occurrence datasets are combined using lines 104 to 159. To access our GBIF downloaded occurrences, use the download codes found in GBIF_codes20210330.txt and the occ_download_get function from the rgbif package.  

-------------------------------

Record conservation status -- "1_Anuran_Conserv_Status.R"

Using the rredlist package, information from the [National Species of Greatest Conservation Need](https://www1.usgs.gov/csas/swap/), and the [federal Endangered Species Act](https://ecos.fws.gov/ecp/services) list, we record the conservation status of the anurans in our species list for those three geopolitical groups. This script produces anuran_conservation_info.csv, which is used as a table key in later scripts. 
- Table S1: Taxonomic and Conservation status for anurans analyzed in this study. Includes number of observations used in analysis
- Table 1: N of species within each conservation status for the three geopolitical groups

-------------------------------

"2_Area of Occupancy Code_Anurans.R"

Code to calculate area of occupancy (AOO) at two grains (small watersheds and 1km Point Buffer) for all species in the final list. Aloso has a for loop to count the observations lost because of the filter by the native range. Output files include  shapefiles at both grains/species, and a CSV with AOO values and ranks for both grain sizes (AOO HUC12 Output_DATE.csv), and an intermediate file that describes the AOO within and outside the contiguous US for both grain sizes (AOOs_raw_DATE.csv). 
- Figure X: Example of occupied watershed across US border and the IUCN native range
- Figure Y: Plot to compare area of watershed polygons within the USGS HUC-12 and HydroBasin datasets

-------------------------------

"3a_ARC_PRISM Download_SCS.R"

This code downloads PRISM data using library(prism). Data may be retrieved through the PRISM website - see https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf for additional details.  

-------------------------------

"3b_Climate Sensitivity_Anurans.R"

This code must be run on the super computer (need higher working memory). Code to calculate standard deviations of climate variables for each species. This requires PRISM rasters downloaded via Script 3a and WorldClim rasters as well. Two large outputs are the average climate values within each occupied polygon (taxa_ws_climate_values_DATE.csv) and six csvs for each species for the climate values within the buffered points (SPP_climate_raster.csv). It also creates the intermediate files that calculate climate sensitivity across the North American range for these species: buff_climate_sums_DATE.csv and range.wide.weighted.CS_DATE.csv.

-------------------------------

"4_RCS_Index_Code.R"

This code calculates Climate Sensitivity (CS) from standard deviation values generated in Script 4. Area of Occupancy (AOO) data (generated in Script 2) and CS data are merged into a single RCS output for all species and grains. Finally, the Relative Climate Sensitivity index (RCS) is calculated for each species and range metric. Three files are created: RCS_Table_DATE.csv, RCS_table_L48_DATE.csv, and the combination of the two former files Anuran RCS values 20210405.csv. 
- Figure S : How AOO varies across spatial extent 
- Figure 4: RCS Dot Plot in ggplot
- Figure S : Sample Size effect on RCS calculation
- Figure S : Distribution of Anuran AOO
- Figure S : AOO dot plot to depict rank vulnerability
- Figure S : Climate Sensitivity dot plot to depict rank vulnerability
- Figure U: Distribution of Climate SDs
- Figure S1: How spatial extent alters RCS calculation - not valid anymore as not all occurrences downloaded (20% cut off)

-------------------------------

"5_RCS_Predict_Conser_Status.R"

Combining the output of 1_Anuran_Conserv_Status and 4_RCS_Index_Code, we calculate differences in mean rank RCS among the different conservation statuses in differen conservation geopolitical groups and spatial extents. We also build a logistic regression (using gams to account for non-linearity to assess if RCS can predict protection status. We also assess how RCS values are distributed among anuran genera
- Figure 2c: Mean rank RCS values and conservation status boxplot
- Figure 2d: RCS poorly predicts protection status 
- Figure 2: RCS values and protection status among Anuran Genera

-------------------------------

"6a_results_paragraph.R"

Code used to summarize and evaluate the results of the above analyses. These summarizations are covered in the results section of the manuscript.

-------------------------------

"6b_Anuran_map_figures.R"

This code aggregates species occurrence at the HUC-12 level to the HUC-8 level. We also calculate the mean RCS value at the HUC-8 level.
- Figure 1a : Anuran species richness at the HUC-8 grain size
- Figure 1b: Average RCS value at the HUC-8 grain size

-------------------------------

"7_anuran_rcs_qaqc.Rmd"

An RMarkdown document recording quality assurance and quality control measures used to evaluate the above analyses. 
- Check that raster::extract correctly computing climate values
- Check that the units match between our two datasets PRISM & WorldClim
- Check for NA values extracted - points that land on the border
- Build a histogram of each species’ climate values (is the SD capturing climate breadth well)



===============================

# Other Files, labelled x:

-------------------------------
RCS_decision_notes.Rmd

This is a document that describes some early decisions when starting this project. It discusses why some taxa are excluded (those with < 20% of their AOO within the contiguous US), why we mostly consider the entire extent, native range extent for the anurans, and a brief comparison of the WorldClim and PRISM datasets. We also explore how the occurrence data shifts through time and how individuals can have marked impacts on certain taxa. 

-------------------------------
make_occ_maps_chk_native_range.Rmd

This creates a document that has the IUCN range maps and occurrence points for each anuran species considered within this analysis. 


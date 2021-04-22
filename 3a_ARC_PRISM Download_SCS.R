# National Fishes Vulnerability Assessment Project - "ARC_PRISM Download_SCS.R"
# Revised by Sam Silknetter, 30June2020

# This code is used to download PRISM climate data for the climate variables of choice. 
# PRISM dataset descriptons: https://prism.oregonstate.edu/documents/PRISM_datasets.pdf

# Install and load necessary libraries.
install.packages('prism')
library(prism)

#Download PRISM Data
options(prism.path = "//home/tracidubose/prism_files/Tmax/")
get_prism_annual(type = "tmax", year = 1895:2019, keepZip = FALSE)

options(prism.path = "/home/tracidubose/prism_files/Tmin/")
get_prism_annual(type = "tmin", year = 1895:2019, keepZip = FALSE)

options(prism.path = "/home/tracidubose/prism_files/ppt/")
get_prism_annual(type = "ppt", year = 1895:2019, keepZip = FALSE)

options(prism.path = "/home/tracidubose/prism_files/Tmean/")
get_prism_annual(type = "tmean", year = 1895:2019, keepZip = FALSE)

options(prism.path = "/home/tracidubose/prism_files/TmaxAug/")
get_prism_monthlys(type = "tmax", year = 1895:2019, mon = 8, keepZip = FALSE)

options(prism.path = "/home/tracidubose/prism_files/TminJan/")
get_prism_monthlys(type = "tmin", year = 1895:2019, mon = 1, keepZip = FALSE)
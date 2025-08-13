#### Download ClimateNA data for GEA analysis


# Install climateNAr package to download data from ClimateNA

#Loading Necessary Packages
library(ClimateNAr)
library(tidyverse)
library(geodata)
#library(raster) #getData will be removed
library(terra)
library(sp)
library(janitor)



#NEW - Note: we migrate to the ClimateNA database due to higher spatial resolution
col_types_list <- cols_only(ID="c", Popln_code="c", Site_code="c", Range_position="f", Latitude="n", Longitude="n", Elevation="n")
input_GEA_coords <- read_csv(file="Rdata/all_popln_geo_coord.csv", col_names = TRUE,  col_select = c(ID, Popln_code, Site_code, Range_position, Longitude, Latitude, Elevation), col_types = col_types_list)


# Derive elevations from WorldClim 30 second (1km^2) database
dem <- rast("Rdata/WorldClim/wc2.1_30s_elev.tif")
points <- vect(input_GEA_coords, geom = c("Longitude", "Latitude"), crs = crs(dem))
input_GEA_coords$Elevation <- extract(dem, points)[, 2] #extract coords
write.csv(input_GEA_coords, file="Rdata/all_popln_geo_coord.csv", row.names=FALSE) #save and replace file



# Download data for GEAs; which uses normals
input_GEA <- input_GEA_coords[, c("Popln_code", "Site_code", "Latitude", "Longitude", "Elevation")]
write.csv(input_GEA, file="Rdata/ClimateNA/input_GEA_coords.csv", row.names=FALSE)
GEA_coords <- "Rdata/ClimateNA/input_GEA_coords.csv"
normal_periods <- c("Normal_1901_1930.nrm", "Normal_1911_1940.nrm", "Normal_1921_1950.nrm",
               "Normal_1931_1960.nrm", "Normal_1941_1970.nrm", "Normal_1951_1980.nrm",
               "Normal_1961_1990.nrm", "Normal_1971_2000.nrm", "Normal_1981_2010.nrm",
               "Normal_1991_2020.nrm")
all_clim_data <- c('YSM') #download yearly, seasonal, and monthly data
outDir="Rdata/ClimateNA/normals/"

climateNAr(inputFile = GEA_coords, periodList = normal_periods, varList = all_clim_data, outDir)

# Download annual data for years 2022-2024, in case this is of interest
yearly_periods <- 2021:2023
all_clim_data <- c('YSM') #download yearly, seasonal, and monthly data
outDir="Rdata/ClimateNA/yearly_2021-2023/"

climateNAr(inputFile = GEA_coords, periodList = yearly_periods, varList = all_clim_data, outDir)
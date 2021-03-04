#############created by Di 2021/02/24##########
####soil properties extracted from HSWD########
rm(list = ls())
setwd("E:/global_soil_database/HWSD")
##update ingestr##
list_pkgs <- c("dplyr", "purrr", "lubridate", "tidyr", "raster", "lubridate", "stringi", "stringr", "sp", "ncdf4", "signal", "climate", "rgdal")
new_pkgs <- list_pkgs[!(list_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)
##update the pachages ingestr needed
#install.packages("lifecycle")
#install.packages("cli")
#install.packages("vctrs")
#install.packages("pillar")
#install.packages("tibble")
#install.packages("cpp11")
#install.packages("BH")
#install.packages("hms")
#install.packages("dplyr")

##update to the latest ingestr
devtools::install_github("stineb/ingestr")

#another way to update ingestr
#if(!require(devtools)){install.packages(devtools)}
#devtools::install_github("stineb/ingestr")
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(purrr)
library(rbeni)
library(ingestr)

## Read data##
df <- read.csv("global_leaf_NP_with_soil_CN_N_water_Al_from_WISE_20210224.csv")%>% 
        rename(lon = lon_estimate,
               lat = lat_estimate,
               elv = Alt_Di_check_final) %>% 
        rename(year_start = Sampling_Year_start,
               year_end = Sampling_Year_end) %>% 
  ## to be done: make sampling year info machine readable: either in the form of "9999" or "8888_9999"
  separate(Sampling_Month, sep = "_", c("month_start", "month_end")) %>% 
  mutate(month_start = as.numeric(month_start), month_end = as.numeric(month_end)) %>% 
  ## arrange by year
  arrange(year_start) %>% 
  ## create identifier
  mutate(id = paste0("i", seq(nrow(.))))%>% 
  ## remove abnormal value [LeafNP>100]
  dplyr::filter(LeafNP<100) 
  ## save with added ids 
write_csv(df, file = "global_leaf_NP_total_Di_20210225_PROC.csv")

## look at some distributions
df %>% 
  ggplot(aes(x = leafN, y = ..density..)) +
  geom_histogram()

df %>% 
  ggplot(aes(x = leafP, y = ..density..)) +
  geom_histogram()

df %>% 
  ggplot(aes(x = LeafNP, y = ..density..)) +
  geom_histogram()

min(df$leafN)
min(df$leafP)
min(df$LeafNP)
max(df$LeafNP)

### Complement missing elevation
df <- read.csv("global_leaf_NP_total_Di_20210225_PROC.csv")

df_tmp <- df %>% 
  dplyr::filter(is.na(elv)) %>% 
  distinct(lon, lat) %>% 
  mutate(sitename = paste0("i", seq(nrow(.))))
##wired: the following ingest is not stable, sometimes need to re-run and then it's ok
df_etopo <- ingest(
  df_tmp,
  source = "etopo1",
  dir = "E:/global_soil_database/HWSD/ETOPO1/"  # adjust this with your local path
) %>% 
  unnest(data) %>% 
  ## add it again so it has lon and lat
  right_join(df_tmp, by = "sitename") %>% 
  ungroup() %>% 
  dplyr::select(elv_etopo = elv, -sitename,lon,lat)

df <- df %>% 
  left_join(df_etopo, by = c("lon", "lat")) %>% 
  rowwise() %>% 
  ungroup()%>% # without this the mutate has an error
  mutate(elv = ifelse(is.na(elv), elv_etopo, elv)) %>% 
  dplyr::select(-elv_etopo) %>% 
  ## create new variable site name (sometimes multiple obs are available for one site)
  mutate(sitename = paste0("i_", as.character(lon), "_", as.character(lat), "_", as.character(elv)))
## save with complemented elevation
write_csv(df, file = "global_leaf_NP_total_Di_20210225_PROC.csv")

## Identify sites: Identify sites (unique lon, lat, elevation)###
df_sites <- df %>%
  distinct(sitename, lon, lat, elv)

## determine start and end year for each site based on available measurements
df_sites <- df %>% dplyr::select(sitename, year_start, year_end) %>% 
  group_by(sitename) %>% 
  summarise(year_start = min(year_start, na.rm = TRUE),
            year_end   = max(year_end, na.rm = TRUE)) %>% 
  right_join(df_sites, by = "sitename")
write_csv(df_sites, file = "df_sites.csv")

## Data overview##
### Observations per species
### How many observations per species?
tmp <- df %>% 
  group_by(Species) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  arrange(desc(count))

tmp %>% 
  ggplot(aes(x = count, y = ..count..)) +
  geom_histogram()

write_csv(tmp, file = "Global_NP_Species_Count.csv")

### Observations per genus###
###How many observations per genus?
tmp <- df %>% 
  group_by(Genus) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  arrange(desc(count))
tmp
tmp %>% 
  ggplot(aes(x = count, y = ..count..)) +
  geom_histogram()
write_csv(tmp, file = "Global_NP_Genus_Count.csv")

### Observations per family
### How many observations per family?
tmp <- df %>% 
  group_by(Family) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  arrange(desc(count))
tmp
tmp %>% 
  ggplot(aes(x = count, y = ..count..)) +
  geom_histogram()
write_csv(tmp, file = "Global_NP_Family_Count.csv")

plot_map_simpl() +
  geom_point(data = df_sites, aes(x = lon, y = lat), color = "red", size = 0.3)

## Identify cells
## bin to half-degree gridcells for determining climate forcing: dlon=0.5, dlat=0.5
## the degree gridcells could be changed here
dlon <- 2
dlat <- 2
lon_breaks <- seq(from = floor(min(df_sites$lon)), to = ceiling(max(df_sites$lon)), by = dlon)
lat_breaks <- seq(from = floor(min(df_sites$lat)), to = ceiling(max(df_sites$lat)), by = dlat)
df_sites <- df_sites %>%
  ungroup() %>% 
  mutate(ilon = cut(lon, 
                    breaks = lon_breaks),
  ilat = cut(lat, 
             breaks = lat_breaks)) %>% 
  mutate(lon_lower = as.numeric( sub("\\((.+),.*", "\\1", ilon)),
         lon_upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", ilon) ),
         lat_lower = as.numeric( sub("\\((.+),.*", "\\1", ilat) ),
         lat_upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", ilat) )) %>% 
  mutate(lon_mid = (lon_lower + lon_upper)/2,
         lat_mid = (lat_lower + lat_upper)/2) %>% 
  ## create cell name to associate with climate input
  mutate(cellname = paste0("icell_", as.character(lon_mid), "_", as.character(lat_mid))) %>% 
  dplyr::select(-ilon, -ilat, -lon_lower, -lon_upper, -lat_lower, -lat_upper)
write_csv(df_sites, file = "df_sites.csv")

df_cells <- df_sites %>% 
  dplyr::select(cellname, lon_mid, lat_mid) %>% 
  distinct()
write_csv(df_cells, file = "df_cells.csv")

## plot sample one site per cell
df_sites_sampled <- df_sites %>% 
  dplyr::select(sitename, lon, lat, elv, cellname) %>%
  group_by(cellname) %>% 
  sample_n(1)

plot_map_simpl() +
  geom_point(data = df_sites_sampled, aes(x = lon, y = lat), color = "red", size = 0.3)

## Complement dataset##
### Climate indeces
##PET must be read separately. It's available from SOFUN outputs. Extract time series and add to `ddf_watch`.
rbeni::extract_pointdata_allsites(
  filename = "global_FULL_MODIS-C006_MOD15A2_v1.4-12-g1b77a9e.a.pet_MEANANN.nc",
  df_lonlat = df_sites %>% slice(1:3))
##This is how it can be done.

## ddf_watch must be a flat data frame (not nested)
## functions to calculate all climate indexes are implemented here:
source("calc_climate_index.R")
## test: use just data for two sites

testcells <- unique(df_sites$sitename)[1:2]

##ddf_watch on stineb should be replaced by df_sites in our code
##there is an error now
df_climate_ind <- df_sites %>% 
  dplyr::filter(sitename %in% testcells) %>% 
  group_by(sitename) %>% 
  nest() %>% 
  mutate(mat = purrr::map(data, ~calc_climate_index_mat(.)),
         matgs = purrr::map(data, ~calc_climate_index_matgs(., temp_base = 5.0)),
         tmonthmin = purrr::map(data, ~calc_climate_index_tmonthmin(.)),
         tmonthmax = purrr::map(data, ~calc_climate_index_tmonthmax(.)),
         ndaysgs = purrr::map(data, ~calc_climate_index_ndaysgs(., temp_base = 5.0)),
         mai = purrr::map(data, ~calc_climate_index_mai(.)),
         maigs = purrr::map(data, ~calc_climate_index_maigs(., temp_base = 5.0)),
         map = purrr::map(data, ~calc_climate_index_map(.)),
         pmonthmin = purrr::map(data, ~calc_climate_index_pmonthmin(.)),
         mapgs = purrr::map(data, ~calc_climate_index_mapgs(., temp_base = 5.0)),
         mavgs = purrr::map(data, ~calc_climate_index_mavgs(., temp_base = 5.0)),
         mav = purrr::map(data, ~calc_climate_index_mav(.))) %>% 
  dplyr::select(-data) %>% 
  unnest(c(mat, matgs, tmonthmin, tmonthmax, ndaysgs, mai, maigs, map, pmonthmin, mapgs, mavgs, mav))




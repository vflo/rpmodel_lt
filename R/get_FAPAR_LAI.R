##### FUNCTIONS FAPAR LAI #######
library(readr)
library(sf)
library(sp)
library(rgdal)
library(tidyverse)
library(data.table)
library(R.utils)
sites_coord <- read_delim("R/data/sites_metadata.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
sites_coord <- sites_coord %>% 
  rename(lon = long,
         lat = lat) %>% 
  dplyr::select(lon, lat)

# Define coordinate reference system
prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my.projection <- st_crs(prj4string)

# Create sf object
lat_long_sf <- st_as_sf(sites_coord, coords = c("lon", "lat"), crs = my.projection)
st_crs(lat_long_sf)

plot(lat_long_sf)

# Export shapefile
st_write(lat_long_sf, "R/data/sites_coord/sites_coord.shp", driver="ESRI Shapefile")


# modification of dist_merge to retain coordinates of y
dist <- function(x1, y1, x2, y2) {
  ((x1-x2)^2 + (y1-y2)^2)^0.5
}

dist_merge_mod <- function(x, y, xeast, xnorth, yeast, ynorth){
  tmp <- t(apply(x[,c(xeast, xnorth)], 1, function(x, y){
    dists <- apply(y, 1, function(x, y) dist(x[2], x[1], y[2], y[1]), x)
    cbind(1:nrow(y), dists)[dists == min(dists),,drop=F][1,]
  }
  , y[,c(yeast, ynorth)]))
  tmp <- cbind(x, 
               min.dist=tmp[,2], 
               y[tmp[,1],-match(c(yeast, ynorth), names(y))],
               y_long = (c(y[tmp[,1], match(c(yeast), names(y))]))[[1]],
               y_lat = (c(y[tmp[,1],match(c(ynorth), names(y))]))[[1]])
  row.names(tmp) <- NULL
  tmp
}

dist_merge <- function(x, y, xeast, xnorth, yeast, ynorth){
  tmp <- t(apply(x[,c(xeast, xnorth)], 1, function(x, y){
    dists <- apply(y, 1, function(x, y) dist(x[2], x[1], y[2], y[1]), x)
    cbind(1:nrow(y), dists)[dists == min(dists),,drop=F][1,]
  }
  , y[,c(yeast, ynorth)]))
  tmp <- cbind(x, min.dist=tmp[,2], y[tmp[,1],-match(c(yeast, ynorth), names(y))])
  row.names(tmp) <- NULL
  tmp
}


sites_coord <- read_delim("R/data/sites_metadata.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

#### MODIS ####
fapar <- read_table("R/data/FAPAR.csv")

fapar %>% 
  mutate(coord = sub('"{""type"":""Point"",""coordinates"":[',"",`.geo`,fixed = TRUE),
         coord = sub(']}"',"",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo

FAPAR_reverse<-dist_merge_mod(sites_coord, foo, 'long', 'lat', 'si_long', 'si_lat')

FAPAR <- foo %>% 
  left_join(FAPAR_reverse%>% 
              dplyr::select(site, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))

write_csv(FAPAR, file="R/data/FAPAR_sites.csv")


#### NOAA ####
fapar <-read_csv("R/data/FAPAR_NOAA.csv")

fapar %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat),
         QC = intToBin(FparLai_QC),
         FparLai_QC = str_sub(QC,-2,-1),
         FparLai_QC = strtoi(FparLai_QC, base = 2))->foo

FAPAR_reverse<-dist_merge_mod(sites_coord, foo, 'long', 'lat', 'si_long', 'si_lat')

FAPAR <- foo %>% 
  dplyr::select(-QC) %>% 
  left_join(FAPAR_reverse%>% 
              dplyr::select(site, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))

write_csv(FAPAR, file="R/data/FAPAR_sites_noaa.csv")

#### ALBEDO ####
albedo <-read_csv("R/data/albedo.csv")

albedo %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo

albedo_reverse<-dist_merge_mod(coord, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

albedo <- foo %>% 
  left_join(albedo_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long)) %>% 
  filter(BRDF_Albedo_Band_Mandatory_Quality_vis == 0)

write_csv(albedo, file="R/data/albedo_sites.csv")
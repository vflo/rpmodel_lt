# install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
#                    "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library(readxl)
library(tidyverse)
library(lmerTest)
library(readr)
library(sf)
library(sp)
library(rgdal)
library(data.table)
library(R.utils)
library("rnaturalearth")
library("rnaturalearthdata")

#### FUNCTIONS ####


calc_dew_t <- function(t,rh){
  if(t<0){
    #Buck, Arden L. 1981. New Equations for Computing Vapor Pressure and Enhancement Factor"
    b = 17.966
    c=247.15
  }else{
    b = 17.368
    c = 238.88
  }
  d = 234.5
  p = exp((b-t/d)*(t/(c+t)))
  c*log(rh/100*p)/(b-log(rh/100*p))
}

calc_rh <- function(t,td){
  if(t<0){
    #Buck, Arden L. 1981. New Equations for Computing Vapor Pressure and Enhancement Factor"
    b = 17.966
    c=247.15
  }else{
    b = 17.368
    c = 238.88
  }
  d = 234.5
  100*exp((b*c*(td-t)/(c+td)+t^2/d)/(c+t))
}

calc_es <- function(t){
  if(t<0){
    #Buck, Arden L. 1981. New Equations for Computing Vapor Pressure and Enhancement Factor"
    b = 17.966
    c=247.15
  }else{
    b = 17.368
    c = 238.88
  }
  a = 611.21 #in Pa
  d = 234.5
  a*exp((b-t/d)*(t/(c+t)))
}

dist <- function(x1, y1, x2, y2) {
  ((x1-x2)^2 + (y1-y2)^2)^0.5
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

# modification of dist_merge to retain coordinates of y
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

#### LEAF DATA LOAD ####
#Wright et al. 2017
leaf_size_df <- read_excel("R/data/aal4760-wright-sm_data_set_s1.xlsx", 
                           sheet = "Global leaf size dataset")

#Schrader et al. 2021
leaf_shape_df <- read_excel("R/data/mcab078_suppl_supplementary_data.xlsx", 
                                               sheet = "Table_S1")%>% 
  separate(species_full, c("genus", "species")) %>% 
  mutate(taxon = paste(genus,species))

#### LEAF WIDTH-AREA RELATIONSHIPS ####

# By Family
mod_width_size_fam <- lm(leaf_area_true~leaf_width:Family,data=leaf_shape_df)
summary(mod_width_size_fam)

# By Genus
mod_width_size_gen<- lm(leaf_area_true~leaf_width:genus,data=leaf_shape_df)
summary(mod_width_size_gen)

# By Species
mod_width_size_sp <- lm(leaf_area_true~leaf_width:taxon,
                         data=leaf_shape_df)
summary(mod_width_size_sp)

# universal
mod_width_size_universal <- lmer(leaf_area_true~leaf_width+(1|Family/genus/taxon),
                        data=leaf_shape_df)
summary(mod_width_size_universal)



#### SITES COORDINATES ####
sites <- leaf_size_df %>% 
  dplyr::select(Site = `Site name`,lat = Latitude, lon = Longitude) %>% 
  group_by(Site) %>% 
  summarise_all(unique)

# # Run only first time
# # Define coordinate reference system
# prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# my.projection <- st_crs(prj4string)
# 
# # Create sf object
# lat_long_sf <- st_as_sf(sites, coords = c("lon", "lat"), crs = my.projection)
# 
# # Checking plot
# theme_set(theme_bw())
# world <- ne_countries(scale = "medium", returnclass = "sf")
# ggplot(data = world) +
#   geom_sf() +
#   geom_point(data=sites,
#              aes(x=lon, y=lat), colour="Deep Pink",
#              fill="Pink",pch=21, size=3, alpha=I(0.7))
# 
# 
# # Save coordinates
# write_sf(lat_long_sf, "R/data/leaf_coord/leaf_size_world_distr.shp")






#### ENVIRONMENTAL VARIABLES DF ####
path_env <- list.files(path="R/data/ERA_values", "*.csv$", full.names=TRUE, recursive = TRUE)
env_df <- purrr::map_df(path_env,read_csv)

    
env_df1 <- env_df %>% 
  rename(
    Tkdew = dew_point,
    Tkair = temperature
    ) %>% 
  rowwise() %>% 
  mutate(
    Tair = Tkair-273.15,
    Tdew = Tkdew-273.15,
    es = calc_es(Tair),
    rh = calc_rh(Tair,Tdew),
    ea = es*rh/100,
    vpd = es-ea,
    ws = sqrt(u_wind_speed^2 + v_wind_speed^2),
    date = lubridate::as_date(Date)
  )


#### FPAR and LAI ####
fapar <- read_csv("R/data/FAPAR_values.csv")%>% 
  dplyr::select(-c(`system:index`,.geo))%>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat),
         QC = intToBin(QC),
         QC = str_sub(QC,-2,-1),
         QC = strtoi(QC, base = 2),
         FAPAR = FAPAR*0.01,
         LAI = LAI*0.1,
         Date = lubridate::ymd(Date)
         )
fapar_coord = fapar %>% dplyr::select(si_lat, si_long) %>% distinct()

fapar_coord <- dist_merge_mod(sites%>% rename(si_lat = lat,si_long = lon), 
                              fapar_coord, 'si_long', 'si_lat', 'si_long', 'si_lat')

fapar <- fapar %>% 
  left_join(fapar_coord %>% 
              dplyr::select(Site,y_lat, y_long) %>%
              rename(si_lat = y_lat,
                     si_long = y_long))%>% 
  filter(QC %in% c(0,2)) %>%
  group_by(Date,Site) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  group_split(Site) %>% 
  purrr::map(function(x){
    z = forecastML::fill_gaps(x[,-5], date_col = 1, 
                              frequency = '1 day',
                              groups = "Site", 
                              static_features = c('si_lat', 'si_long'))
    z <- z %>% mutate(FAPAR = as.numeric(FAPAR),
                      LAI = as.numeric(LAI))
    FAPAR = zoo::na.approx(z[,c(3)],na.rm=FALSE, maxgap = as.numeric(100))
    LAI = zoo::na.approx(z[,c(4)],na.rm=FALSE, maxgap = as.numeric(100))
    z$FAPAR <- FAPAR
    z$LAI <- LAI

    z
  }) %>% bind_rows()


# env_df <- env_df %>% left_join(fapar)

emi <- read_csv("R/data/EMISSIVITY_values.csv")%>% 
  dplyr::select(-c(`system:index`,.geo)) %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat),
         Emis_31 = 0.49+Emis_31*0.002,
         Emis_32 = 0.49+Emis_32*0.002,
         QC_day = intToBin(QC_day),
         QC_day = str_sub(QC_day,-2,-1),
         QC_day = strtoi(QC_day, base = 2),
         QC_night = intToBin(QC_night),
         QC_night = str_sub(QC_night,-2,-1),
         QC_night = strtoi(QC_night, base = 2),
         Date = lubridate::ymd(Date)
         )
emi_coord <- emi %>% dplyr::select(si_lat, si_long) %>% distinct()

emi_coord <- dist_merge_mod(sites%>% rename(si_lat = lat,si_long = lon), 
                            emi_coord, 'si_long', 'si_lat', 'si_long', 'si_lat')

emi <- emi %>% 
  left_join(emi_coord %>% 
              dplyr::select(Site,y_lat, y_long) %>%
              rename(si_lat = y_lat,
                     si_long = y_long))%>% 
  filter(QC_day %in% c(0)) %>%
  group_by(Date,Site) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  group_split(Site) %>% 
  purrr::map(function(x){
    z = forecastML::fill_gaps(x[,-5], date_col = 1, 
                              frequency = '1 day',
                              groups = "Site", 
                              static_features = c('si_lat', 'si_long'))
    z <- z %>% mutate(Emis_31 = as.numeric(Emis_31),
                      Emis_32 = as.numeric(Emis_32))
    Emis_31 = zoo::na.approx(z[,c(3)],na.rm=FALSE, maxgap = as.numeric(100))
    Emis_32 = zoo::na.approx(z[,c(4)],na.rm=FALSE, maxgap = as.numeric(100))
    z$Emis_31 <- Emis_31
    z$Emis_32 <- Emis_32

    z
  }) %>% bind_rows()

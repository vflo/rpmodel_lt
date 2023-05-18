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

# This function calculates the dew point temperature from air temperature and relative humidity
# using the Buck- Arden L. equation for vapor pressure and enhancement factor.
calc_dew_t <- function(t, rh){
  ifelse(t<0, {
    b = 17.966
    c = 247.15
  }, {
    b = 17.368
    c = 238.88
  })
  d = 234.5
  p = exp((b-t/d)*(t/(c+t)))
  c*log(rh/100*p)/(b-log(rh/100*p))
}

# This function calculates the relative humidity from air temperature and dew point temperature
# using the Buck- Arden L. equation.
calc_rh <- function(t, td) {
  b <- ifelse(t < 0, 17.966, 17.368)
  c <- ifelse(t < 0, 247.15, 238.88)
  d <- 234.5
  100 * exp((b * c * (td - t) / (c + td) + t^2 / d) / (c + t))
}

# This function calculates the saturation vapor pressure from air temperature
# using the Buck- Arden L. equation.
calc_es <- function(t_vec) {
  b_vec <- ifelse(t_vec < 0, 17.966, 17.368)
  c_vec <- ifelse(t_vec < 0, 247.15, 238.88)
  a <- 611.21 # in Pa
  d <- 234.5
  a * exp((b_vec - t_vec / d) * (t_vec / (c_vec + t_vec)))
}

# This function calculates the distance between two points in a 2D coordinate system.
dist <- function(x1, y1, x2, y2) {
  ((x1-x2)^2 + (y1-y2)^2)^0.5
}

# This function takes two data frames with x and y coordinates and merges them, calculating the
# closest point in the second data frame to each point in the first data frame and adding this
# information to a new column of the first data frame.
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
  dplyr::group_by(Site) %>% 
  dplyr::summarise(across(everything(),unique))
time_zone = lutz::tz_lookup_coords(sites$lat, sites$lon,"accurate")
timezones <- tibble(Site = sites$Site, time_zone = time_zone)

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
env_df <- purrr::map_df(path_env,read_csv)%>% 
  dplyr::select(-c(`system:index`,.geo, Hour)) %>% 
  na.omit()

env_df_coord = env_df %>% dplyr::select(si_lat, si_long) %>% distinct()

env_df_coord <- dist_merge_mod(sites%>% rename(si_lat = lat,si_long = lon), 
                              env_df_coord, 'si_long', 'si_lat', 'si_long', 'si_lat')

    
env_df <- env_df_coord %>%
  dplyr::select(Site, y_lat, y_long) %>%
  dplyr::rename(si_lat = y_lat, si_long = y_long) %>%
  dplyr::left_join(env_df, relationship = "many-to-many") %>% 
  dplyr::left_join(timezones)

env_df$Tair = env_df$temperature - 273.15
env_df$Tdew = env_df$dew_point - 273.15
env_df$es = calc_es(env_df$Tair)
env_df$rh = calc_rh(env_df$Tair, env_df$Tdew)
env_df$ea = env_df$es * env_df$rh / 100
env_df$vpd = env_df$es - env_df$ea
env_df$ws = sqrt(env_df$u_wind_speed ^ 2 + env_df$v_wind_speed ^ 2)
env_df$sw_in = env_df$sw_in / 3600
env_df$sw_in_net = env_df$sw_in_net / 3600
env_df$ppfd = env_df$sw_in * 4.6 * 0.5
env_df$timestamp = lubridate::ymd_hms(env_df$Date)
# env_df$timestamp =  purrr::map(env_df %>% split(seq(nrow(.))),
#                                function(x) {
#                                  lubridate::with_tz(x$timestamp_origin,x$time_zone)
#                                  }) %>% bind_rows()
env_df$Date = lubridate::as_date(env_df$timestamp)
env_df <- env_df%>%
  dplyr::distinct()



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
    FAPAR = ifelse(FAPAR<0,0,FAPAR)
    FAPAR = ifelse(FAPAR>1,1,FAPAR)
    LAI = ifelse(LAI<0,0,LAI)
    z$FAPAR <- FAPAR
    z$LAI <- LAI

    z
  }) %>% bind_rows()




#### EMISSIVITY ####
emi <- read_csv("R/data/EMISSIVITY_values.csv") %>% 
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
                     si_long = y_long),
            relationship = "many-to-many") %>% 
  filter(QC_day %in% c(0,1)) %>%
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
    
    z[is.na(z$Emis_31),"Emis_31"] <- mean(z$Emis_31,na.rm = TRUE)
    z[is.na(z$Emis_32),"Emis_32"] <- mean(z$Emis_32,na.rm = TRUE)


    z
  }) %>% bind_rows()




#### VEGETATION HEIGHT ####
veg_height <- read_csv("R/data/veg_height.csv") 

veg_height_coord <- veg_height %>% dplyr::select(si_lat, si_long) %>% distinct()

veg_height_coord <- dist_merge_mod(sites%>% rename(si_lat = lat,si_long = lon), 
                            veg_height_coord, 'si_long', 'si_lat', 'si_long', 'si_lat')

veg_height <- veg_height %>% 
  left_join(veg_height_coord %>% 
              dplyr::select(Site,y_lat, y_long) %>%
              rename(si_lat = y_lat,
                     si_long = y_long)) %>% 
  distinct() %>% 
  mutate(veg_height = case_when(is.na(mean)~1,
                                TRUE~`mean`)) %>% 
  dplyr::select(-mean)




#### CO2 ####
co2 <- read_csv("R/data/co2_mm_mlo.csv", 
                   skip = 56)  %>% 
  mutate(Date = lubridate::ym(paste0(year,"-",month)))%>% 
  dplyr::select(Date,average) %>% 
  dplyr::rename(co2 = average)
z = forecastML::fill_gaps(co2, date_col = 1, 
                          frequency = '1 day',
                          groups = NULL, 
                          static_features = NULL)
co2 <- zoo::na.approx(z[,2])
z$co2 <- co2
co2 <- z

rm(z)

#### ELEVATION ####
sites_elevation <- leaf_size_df %>% 
  dplyr::select(Site = `Site name`, elev = `Elevation (m)`) %>%
  group_by(Site) %>% 
  summarise(across(everything(),unique))


#### FINAL DATASET ####
env_df <- env_df %>% 
  # filter(Site %in% c("Abisko","Abrams_Pennsylvania")) %>%
  left_join(fapar %>% dplyr::select(-c(si_lat,si_long)) , by=c("Site","Date")) %>% 
  left_join(emi %>% dplyr::select(-c(si_lat,si_long)), by=c("Site","Date")) %>% 
  left_join(co2, by=c("Date")) %>% 
  left_join(veg_height %>% dplyr::select(-c(si_lat,si_long)), by=c("Site")) %>% 
  left_join(sites_elevation, by=c("Site"))

env_df[is.na(env_df$Emis_31),"Emis_31"] <- 0.97
env_df[is.na(env_df$Emis_32),"Emis_32"] <- 0.97


env_df %>% 
  group_by(Site) %>% 
  group_split() %>% 
  purrr::map(function(x){
    write_csv(x,paste0("R/data/df_opt/", gsub("/","_",unique(x$Site)),".csv"))
  })




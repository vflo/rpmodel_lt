#### TEST SCRIPT ####

########################################################################
# 01.load the libraries
########################################################################
library(xts)
# library(Evapotranspiration)
# library(readFlux)
library(bigleaf)
library(ggplot2)
library(zoo)
library(tidyverse)
if(!require(devtools)){install.packages("devtools")}
# devtools::install_github("aukkola/FluxnetLSM", build_vignettes = TRUE)
# library("FluxnetLSM")
sapply(list("R/rpmodel_core.R","R/rpmodel.R","R/rpmodel_subdaily.R",
            "R/subroutines.R","R/include_fapar_lai.R","R/include_albedo.R"),source,.GlobalEnv)
library(MASS)
<<<<<<< HEAD
library(ggfx)
=======
library(ggplot2)
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
library(viridis)
library(Metrics)
library(sf)
library(R.utils)
library(rsplashtest)
source("R/evaluation_sim_plot.R")
<<<<<<< HEAD

library(lme4)
library(lmerTest)
library(relaimpo)
library(DAAG)
library(party)
library(rpart)
library(rpart.plot)
library(mlbench)
library(caret)
library(pROC)
library(tree)
library(ggplot2)
library(ggpubr)
library(effects)
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)
library(tmap)
library(readxl)
library(ggfx)
library(caret)
source('R/sfn_predict_nn.R')

create_train_test <- function(data, size = 0.8) {
  ind <- sample(2, nrow(data), replace = T, prob = c(size, 1-size))
  return (list(data[ind == 1, ],data[ind == 2, ]))
}

accuracy_tune <- function(fit) {
  predict_unseen <- predict(fit, data_test, type = 'class')
  table_mat <- table(data_test$cat_zero, predict_unseen)
  accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
  accuracy_Test
}

=======
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Summarise bbox of soil grids
summary_soil_bbox <- function(foo, aggregation){
  soil <- data_frame(widths = c(foo[[1]]$widths),
                     clay = c(mean(c(foo[[1]][1,2],
                                     foo[[2]][1,2],
                                     foo[[3]][1,2],
                                     foo[[4]][1,2],
                                     foo[[5]][1,2]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][2,2],
                                     foo[[2]][2,2],
                                     foo[[3]][2,2],
                                     foo[[4]][2,2],
                                     foo[[5]][2,2]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][3,2],
                                     foo[[2]][3,2],
                                     foo[[3]][3,2],
                                     foo[[4]][3,2],
                                     foo[[5]][3,2]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][4,2],
                                     foo[[2]][4,2],
                                     foo[[3]][4,2],
                                     foo[[4]][4,2],
                                     foo[[5]][4,2]),
                                   na.rm  = TRUE)),
                     sand = c(mean(c(foo[[1]][1,3],
                                     foo[[2]][1,3],
                                     foo[[3]][1,3],
                                     foo[[4]][1,3],
                                     foo[[5]][1,3]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][2,3],
                                     foo[[2]][2,3],
                                     foo[[3]][2,3],
                                     foo[[4]][2,3],
                                     foo[[5]][2,3]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][3,3],
                                     foo[[2]][3,3],
                                     foo[[3]][3,3],
                                     foo[[4]][3,3],
                                     foo[[5]][3,3]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][4,3],
                                     foo[[2]][4,3],
                                     foo[[3]][4,3],
                                     foo[[4]][4,3],
                                     foo[[5]][4,3]),
                                   na.rm  = TRUE)),
                     om =   c(mean(c(foo[[1]][1,4],
                                     foo[[2]][1,4],
                                     foo[[3]][1,4],
                                     foo[[4]][1,4],
                                     foo[[5]][1,4]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][2,4],
                                     foo[[2]][2,4],
                                     foo[[3]][2,4],
                                     foo[[4]][2,4],
                                     foo[[5]][2,4]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][3,4],
                                     foo[[2]][3,4],
                                     foo[[3]][3,4],
                                     foo[[4]][3,4],
                                     foo[[5]][3,4]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][4,4],
                                     foo[[2]][4,4],
                                     foo[[3]][4,4],
                                     foo[[4]][4,4],
                                     foo[[5]][4,4]),
                                   na.rm  = TRUE)),
                     bd =   c(mean(c(foo[[1]][1,5],
                                     foo[[2]][1,5],
                                     foo[[3]][1,5],
                                     foo[[4]][1,5],
                                     foo[[5]][1,5]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][2,5],
                                     foo[[2]][2,5],
                                     foo[[3]][2,5],
                                     foo[[4]][2,5],
                                     foo[[5]][2,5]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][3,5],
                                     foo[[2]][3,5],
                                     foo[[3]][3,5],
                                     foo[[4]][3,5],
                                     foo[[5]][3,5]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][4,5],
                                     foo[[2]][4,5],
                                     foo[[3]][4,5],
                                     foo[[4]][4,5],
                                     foo[[5]][4,5]),
                                   na.rm  = TRUE)),
                     rfc =  c(mean(c(foo[[1]][1,6],
                                     foo[[2]][1,6],
                                     foo[[3]][1,6],
                                     foo[[4]][1,6],
                                     foo[[5]][1,6]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][2,6],
                                     foo[[2]][2,6],
                                     foo[[3]][2,6],
                                     foo[[4]][2,6],
                                     foo[[5]][2,6]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][3,6],
                                     foo[[2]][3,6],
                                     foo[[3]][3,6],
                                     foo[[4]][3,6],
                                     foo[[5]][3,6]),
                                   na.rm  = TRUE),
                              mean(c(foo[[1]][4,6],
                                     foo[[2]][4,6],
                                     foo[[3]][4,6],
                                     foo[[4]][4,6],
                                     foo[[5]][4,6]),
                                   na.rm  = TRUE)),
                     aggregation = aggregation
  )
  return(soil)
  
}


get_soil <- function(coord_object){
  soil <- medfateutils::soilgridsParams(coord_object,
                                        widths = 1000)[1,]
  soil$aggregation <- 0
  if(any(is.na(soil))){
    foo <- medfateutils::soilgridsParams(coord_object %>%
                                           st_buffer(250) %>%
                                           st_bbox() %>%
                                           st_as_sfc(),
                                         widths = 1000
    )
    soil <- summary_soil_bbox(foo,250)[1,]
  }
  if(any(is.na(soil))){
    foo <- medfateutils::soilgridsParams(coord_object %>%
                                           st_buffer(500) %>%
                                           st_bbox() %>%
                                           st_as_sfc(),
                                         widths = 1000
    )
    soil <- summary_soil_bbox(foo,500)[1,]
  }
  if(any(is.na(soil))){
    foo <- medfateutils::soilgridsParams(coord_object %>%
                                           st_buffer(1000) %>%
                                           st_bbox() %>%
                                           st_as_sfc(),
                                         widths = 1000
    )
    soil <- summary_soil_bbox(foo,1000)[1,]
  }
  if(any(is.na(soil))){
    foo <- medfateutils::soilgridsParams(coord_object %>%
                                           st_buffer(4000) %>%
                                           st_bbox() %>%
                                           st_as_sfc(),
                                         widths = 1000
    )
    soil <- summary_soil_bbox(foo,5000)[1,]
  }
  if(any(is.na(soil))){
    foo <- medfateutils::soilgridsParams(coord_object %>%
                                           st_buffer(10000) %>%
                                           st_bbox() %>%
                                           st_as_sfc(),
                                         widths = 1000
    )
    soil <- summary_soil_bbox(foo,10000)[1,]
  }
  soil$site <- coord_object$site
  return(soil)
}



### splash wrapper
runsplash<-function(df,lat,elev,slop,asp,au,soil){
  #### whc of 150mm ~ deph<1m
  soil_data <- soil#c(unique(soil$sand),unique(soil$clay),unique(soil$om),unique(soil$rfc),unique(soil$bd),1)
  prec <- df %>% group_by(Date) %>% dplyr::select(total_precipitation_hourly) %>% 
    summarise(prec = sum(total_precipitation_hourly,na.rm=TRUE))
  sw_in <- df %>% group_by(Date) %>% dplyr::select(sw_in) %>% summarise(sw_in = mean(sw_in,na.rm=TRUE))
  ta <- df %>% group_by(Date) %>% dplyr::select(Tair) %>% summarise(ta = mean(Tair,na.rm=TRUE))
  as.data.frame(
    splash.point(
      sw_in=sw_in$sw_in,	# shortwave radiation W/m2
      tc=ta$ta,		# air temperature C
      pn= prec$prec*1000,		# precipitation mm
      lat=lat,		# latitude deg
      elev=elev,		# elevation masl
      slop=slop,	# slope deg
      asp=asp,		# aspect deg
      soil_data=soil_data, 		# soil data: sand,clay,som in w/w %. Gravel v/v %, bulk density g/cm3, and depth to the bedrock (m)**
      Au=au,		# upslope area m2
      resolution=90.0, # resolution pixel dem used to get Au
      time_index=ta$Date,
      ts_out=FALSE
    )
  )
}



#### DATA LOAD ####
sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")




<<<<<<< HEAD
precip_ERA5 <- read_csv(path="R/data/ERA_values_prec_amf.csv")%>% 
=======
precip_ERA5 <- read_csv(path="R/data/ERA_values_prec_amf.csv$")%>% 
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
  dplyr::select(-c(`system:index`,.geo, Hour)) %>% 
  na.omit()

env_df_coord = precip_ERA5 %>% dplyr::select(si_lat, si_long) %>% distinct()

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

 
# 
# lat_long_sf <- sf::st_as_sf(data.table::data.table(sites_metadata[c(12,14,15,31),]),
#                             coords = c("long","lat"), crs = 4326)
# 
# sf::st_write(lat_long_sf, "R/data/amf_coord/amf_coord.shp")
# soil <- list()
# for(i in 1:nrow(lat_long_sf)
# ){
#   soil[[i]] <- get_soil(lat_long_sf[i,])
#   
# }
# soil <- bind_rows(soil)
# 
# save(soil, file="R/data/soil_amf.Rdata")
load(file="R/data/soil_amf.Rdata")

#### EMISSIVITY ####
emi <- read_csv("R/data/EMISSIVITY_values_AMF.csv") %>% 
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


# 01.find the file name ------------------------------------------------

#load files available
#### read the files' paths
filenames.fluxnet<- list.files(path="R/data/final_sites", "*.csv$", full.names=TRUE,recursive = TRUE)


# 02.Read the data US-Ha -----------------------------------------------

filename <- filenames.fluxnet[10]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

data_flx_HA <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(P) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP),
         P = P/1000)
summary(data_flx_HA)

env_df <- data_flx_HA %>% 
  # filter(Site %in% c("Abisko","Abrams_Pennsylvania")) %>%
  left_join(soil , by=c("site")) %>% 
  rename(Date = date,
         total_precipitation_hourly = P,
         sw_in = SW_IN)

foo <- runsplash(env_df,unique(env_df$lat),unique(env_df$elev), 
                 slop = 0, asp = 0, au = 0, 
                 soil =  c(unique(env_df$sand),unique(env_df$clay),
                           unique(env_df$om),unique(env_df$rfc),
                           unique(env_df$bd),1))
foo <- foo %>% bind_cols( env_df %>% dplyr::select(Date) %>% distinct())
data_flx_HA <- left_join(env_df,foo)
AI = data_flx_HA %>% 
  mutate(year = lubridate::year(Date)) %>%
  group_by(year) %>% 
  summarise(prec = sum(total_precipitation_hourly,na.rm=TRUE),
            pet = sum(pet,na.rm=TRUE)/24/1000) %>% 
  ungroup() %>% 
  summarise(prec = mean(prec,na.rm=TRUE),
            pet = mean(pet,na.rm=TRUE),
            AI = pet/prec)

gc()

data_flx_HA <- data_flx_HA  %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar)


res_HA <- rpmodel_subdaily(TIMESTAMP = data_flx_HA$timestamp, tc = data_flx_HA$Tair, vpd = data_flx_HA$VPD*100,
                           co2 = data_flx_HA$CO2, fapar = data_flx_HA$FAPAR, LAI = data_flx_HA$LAI,
                           ppfd = data_flx_HA$PPFD, u=data_flx_HA$WS, ustar = NA,# data_flx_HA$Ustar, 
                           canopy_height = data_flx_HA$Hc, patm =data_flx_HA$PA*1000, 
                           elv = data_flx_HA$elev, z = data_flx_HA$z, leafwidth = 0.05, 
                           netrad = NA, #data_flx_HA$NETRAD,
                           beta = 146.0, c_cost = 0.41, soilm = data_flx_HA$sm_lim, AI = AI$AI,
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, do_soilmstress = TRUE,
                           upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE, epsleaf = 0.98, #thermal absorptivity of the leaf
                           energy_params = list(
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_HA_swc <- as_tibble(res_HA) %>% cbind(data_flx_HA)

save(df_res_HA_swc,file="R/data/df_res_HA_swc.RData")
load("R/data/df_res_HA_swc.RData")
load("R/data/df_res_HA.RData")
df_res_HA <- df_res_HA_swc

df_res_HA %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  pivot_longer(c(d_can_air,d_can_leaf)) %>%
  ggplot() +
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_smooth(mapping = aes(hour,value, color=name, fill=name),linewidth = 1, 
              show.legend = FALSE)+
  geom_point(aes(hour,value, color=name),alpha=0.20)+
  facet_wrap(~month)+
  theme_bw()+xlab("Hour")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                'd_can_leaf' = expression("("*T[sim] - T[air]*")")),
                     values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top")+
  NULL

df_res_HA %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*100,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("VPD [Pa]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_HA %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(Tair,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("T air [ºC]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_HA %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric(),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(WS,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("WS [m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_HA %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(PPFD,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("PPFD ["*mu*"mol"*m^{-2}~s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL


df_res_HA %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*WS,d_can_air,color=PPFD))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("VPD x WS [Pa m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL


HA_1 <- df_res_HA %>%
  filter(PPFD>0) %>% 
  ggplot(aes(lE,LE))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 500, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_HA$lE, df_res_HA$LE,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

HA_2 <- df_res_HA %>%
  ggplot(aes(Rnet, NETRAD))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 662, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_HA$Rnet, df_res_HA$NETRAD,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

HA_3 <- df_res_HA %>%
  ggplot(aes(tcleaf, Tcan))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = -10, y = 23.7, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_HA$tcleaf, df_res_HA$Tcan,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Canopy temperature [ºC]"))+
  xlab(expression("Simulated Leaf temperature [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

ggpubr::ggarrange(HA_1,HA_2,HA_3, labels = c("a","b","c"),ncol = 3)
#T
# plot.eval.dens(df_res_HA$tcleaf,df_res_HA$Tcan)





# 03.Read the data US-Me2 ----------------------------------------------

filename <- filenames.fluxnet[12]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

data_flx_ME <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(P,Tair,SW_IN) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP),
         P = P/1000)
summary(data_flx_ME)

env_df <- data_flx_ME %>% 
  # filter(Site %in% c("Abisko","Abrams_Pennsylvania")) %>%
  left_join(soil , by=c("site")) %>% 
  rename(Date = date,
         total_precipitation_hourly = P,
         sw_in = SW_IN)

foo <- runsplash(env_df, unique(env_df$lat),unique(env_df$elev), 
                 slop = 0, asp = 0, au = 0, 
                 soil =  c(unique(env_df$sand),unique(env_df$clay),
                           unique(env_df$om),unique(env_df$rfc),
                           unique(env_df$bd),1))
foo <- foo %>% bind_cols( env_df %>% dplyr::select(Date) %>% distinct())
data_flx_ME <- left_join(env_df,foo)
AI = data_flx_ME %>% 
  mutate(year = lubridate::year(Date)) %>%
  group_by(year) %>% 
  summarise(prec = sum(total_precipitation_hourly,na.rm=TRUE),
            pet = sum(pet,na.rm=TRUE)/24/1000) %>% 
  ungroup() %>% 
  summarise(prec = mean(prec,na.rm=TRUE),
            pet = mean(pet,na.rm=TRUE),
            AI = pet/prec)

gc()

data_flx_ME <- data_flx_ME  %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar)

res_ME <- rpmodel_subdaily(TIMESTAMP = data_flx_ME$timestamp, tc = data_flx_ME$Tair, vpd = data_flx_ME$VPD*100,
                           co2 = data_flx_ME$CO2, fapar = data_flx_ME$FAPAR, LAI = data_flx_ME$LAI,
                           ppfd = data_flx_ME$PPFD, u=data_flx_ME$WS, ustar = NA,# data_flx_ME$Ustar, 
                           canopy_height = data_flx_ME$Hc, patm =data_flx_ME$PA*1000, 
                           elv = data_flx_ME$elev, z = data_flx_ME$z, leafwidth = 0.001, 
                           netrad = NA, #data_flx_ME$NETRAD,
                           beta = 146.0, c_cost = 0.41, soilm = data_flx_ME$sm_lim, AI = AI$AI,
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, do_soilmstress = TRUE,
                           upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE, epsleaf = 0.98, #thermal absorptivity of the leaf
                           energy_params = list(
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_ME_swc <- as_tibble(res_ME) %>% cbind(data_flx_ME)

save(df_res_ME_swc,file="R/data/df_res_ME_swc.RData")
load("R/data/df_res_ME.RData")
df_res_ME <- df_res_ME_swc

df_res_ME %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  pivot_longer(c(d_can_air,d_can_leaf)) %>%
  ggplot() +
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_smooth(mapping = aes(hour,value, color=name, fill=name),linewidth = 1, 
              show.legend = FALSE)+
  geom_point(aes(hour,value, color=name),alpha=0.20)+
  facet_wrap(~month)+
  theme_bw()+xlab("Hour")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                'd_can_leaf' = expression("("*T[sim] - T[air]*")")),
                     values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top")+
  NULL

df_res_ME %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*100,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("VPD [Pa]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_ME %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(Tair,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("T air [ºC]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_ME %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric(),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(WS,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("WS [m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_ME %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(PPFD,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("PPFD ["*mu*"mol"*m^{-2}~s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL


df_res_ME %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*WS,d_can_air,color=PPFD))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("VPD x WS [Pa m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL


ME_1 <- df_res_ME %>%
  filter(PPFD>0) %>% 
  ggplot(aes(lE,LE))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 500, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_ME$lE, df_res_ME$LE,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

ME_2 <- df_res_ME %>%
  ggplot(aes(Rnet, NETRAD))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 760, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_ME$Rnet, df_res_ME$NETRAD,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

ME_3 <- df_res_ME %>%
  ggplot(aes(tcleaf, Tcan))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 0, y = 28.8, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_ME$tcleaf, df_res_ME$Tcan,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Canopy temperature [ºC]"))+
  xlab(expression("Simulated Leaf temperature [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

ggpubr::ggarrange(ME_1,ME_2,ME_3, labels = c("a","b","c"),ncol = 3)
#T
# plot.eval.dens(df_res_ME$tcleaf,df_res_ME$Tcan)

plot.eval.dens(df_res_ME$lE,df_res_ME$LE)


# 04.Read the data US-NR1 ----------------------------------------------

filename <- filenames.fluxnet[13]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

data_flx_NR <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(P,Tair,SW_IN) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP),
         P = P/1000)
summary(data_flx_NR)

env_df <- data_flx_NR %>% 
  # filter(Site %in% c("Abisko","Abrams_Pennsylvania")) %>%
  left_join(soil , by=c("site")) %>% 
  rename(Date = date,
         total_precipitation_hourly = P,
         sw_in = SW_IN)

foo <- runsplash(env_df,unique(env_df$lat),unique(env_df$elev), 
                 slop = 0, asp = 0, au = 0, 
                 soil =  c(unique(env_df$sand),unique(env_df$clay),
                           unique(env_df$om),unique(env_df$rfc),
                           unique(env_df$bd),1))
foo <- foo %>% bind_cols( env_df %>% dplyr::select(Date) %>% distinct())
data_flx_NR <- left_join(env_df,foo)
AI = data_flx_NR %>% 
  mutate(year = lubridate::year(Date)) %>%
  group_by(year) %>% 
  summarise(prec = sum(total_precipitation_hourly,na.rm=TRUE),
            pet = sum(pet,na.rm=TRUE)/24/1000) %>% 
  ungroup() %>% 
  summarise(prec = mean(prec,na.rm=TRUE),
            pet = mean(pet,na.rm=TRUE),
            AI = pet/prec)

gc()

data_flx_NR <- data_flx_NR  %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar)

res_NR <- rpmodel_subdaily(TIMESTAMP = data_flx_NR$timestamp, tc = data_flx_NR$Tair, vpd = data_flx_NR$VPD*100,
                           co2 = data_flx_NR$CO2, fapar = data_flx_NR$FAPAR, LAI = data_flx_NR$LAI,
                           ppfd = data_flx_NR$PPFD, u=data_flx_NR$WS, ustar = NA,# data_flx_NR$Ustar, 
                           canopy_height = data_flx_NR$Hc, patm =data_flx_NR$PA*1000, 
                           elv = data_flx_NR$elev, z = data_flx_NR$z, leafwidth = 0.001, 
                           netrad = NA, #data_flx_NR$NETRAD,
                           beta = 146.0, c_cost = 0.41, soilm = data_flx_NR$sm_lim, AI = AI$AI,
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, do_soilmstress = TRUE,
                           upscaling_method = "daily", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE, epsleaf = 0.98, #thermal absorptivity of the leaf
                           energy_params = list(
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_NR <- as_tibble(res_NR) %>% cbind(data_flx_NR)

save(df_res_NR,file="R/data/df_res_NR.RData")

df_res_NR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  pivot_longer(c(d_can_air,d_can_leaf)) %>%
  ggplot() +
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_smooth(mapping = aes(hour,value, color=name, fill=name),linewidth = 1, 
              show.legend = FALSE)+
  geom_point(aes(hour,value, color=name),alpha=0.20)+
  facet_wrap(~month)+
  theme_bw()+xlab("Hour")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                'd_can_leaf' = expression("("*T[sim] - T[air]*")")),
                     values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top")+
  NULL

df_res_NR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*100,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("VPD [Pa]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_NR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(Tair,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("T air [ºC]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_NR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric(),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(WS,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("WS [m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_NR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(PPFD,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("PPFD ["*mu*"mol"*m^{-2}~s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL


df_res_NR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*WS,d_can_air,color=PPFD))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("VPD x WS [Pa m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

NR_1 <- df_res_NR %>%
  filter(PPFD>0,LE>0) %>% 
  ggplot(aes(lE,LE))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 182, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_NR$lE, df_res_NR$LE,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

NR_2 <- df_res_NR %>%
  ggplot(aes(Rnet, NETRAD))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 751, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_NR$Rnet, df_res_NR$NETRAD,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

NR_3 <- df_res_NR %>%
  ggplot(aes(tcleaf, Tcan))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = -10, y = 16.3, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_NR$tcleaf, df_res_NR$Tcan,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Canopy temperature [ºC]"))+
  xlab(expression("Simulated Leaf temperature [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

ggpubr::ggarrange(NR_1,NR_2,NR_3, labels = c("a","b","c"),ncol = 3)
#T
plot.eval.dens(df_res_NR$tcleaf,df_res_NR$Tcan)





# 05.Read the data US-Wrc ----------------------------------------------

filename <- filenames.fluxnet[29]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)


data_flx_WR <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(P,Tair,SW_IN) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP),
         P = P/1000)
summary(data_flx_WR)

env_df <- data_flx_WR %>% 
  # filter(Site %in% c("Abisko","Abrams_Pennsylvania")) %>%
  left_join(soil , by=c("site")) %>% 
  rename(Date = date,
         total_precipitation_hourly = P,
         sw_in = SW_IN)

foo <- runsplash(env_df,unique(env_df$lat),unique(env_df$elev), 
                 slop = 0, asp = 0, au = 0, 
                 soil =  c(unique(env_df$sand),unique(env_df$clay),
                           unique(env_df$om),unique(env_df$rfc),
                           unique(env_df$bd),1))
foo <- foo %>% bind_cols( env_df %>% dplyr::select(Date) %>% distinct())
data_flx_WR <- left_join(env_df,foo)
AI = data_flx_WR %>% 
  mutate(year = lubridate::year(Date)) %>%
  group_by(year) %>% 
  summarise(prec = sum(total_precipitation_hourly,na.rm=TRUE),
            pet = sum(pet,na.rm=TRUE)/24/1000) %>% 
  ungroup() %>% 
  summarise(prec = mean(prec,na.rm=TRUE),
            pet = mean(pet,na.rm=TRUE),
            AI = pet/prec)

gc()

data_flx_WR <- data_flx_WR  %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar)

res_WR <- rpmodel_subdaily(TIMESTAMP = data_flx_WR$timestamp, tc = data_flx_WR$Tair, vpd = data_flx_WR$VPD*100,
                           co2 = data_flx_WR$CO2, fapar = data_flx_WR$FAPAR, LAI = data_flx_WR$LAI,
                           ppfd = data_flx_WR$PPFD, u=data_flx_WR$WS, ustar = NA,# data_flx_WR$Ustar, 
                           canopy_height = data_flx_WR$Hc, patm =data_flx_WR$PA*1000, 
                           elv = data_flx_WR$elev, z = data_flx_WR$z, leafwidth = 0.001, 
                           netrad = NA, #data_flx_WR$NETRAD,
                           beta = 146.0, c_cost = 0.41, soilm = data_flx_WR$sm_lim, AI = AI$AI,
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, do_soilmstress = TRUE,
                           upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE,epsleaf = 0.98,
                           energy_params = list(
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_WR <- as_tibble(res_WR) %>% cbind(data_flx_WR)

save(df_res_WR,file="R/data/df_res_WR.RData")

df_res_WR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  pivot_longer(c(d_can_air,d_can_leaf)) %>%
  ggplot() +
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_smooth(mapping = aes(hour,value, color=name, fill=name),linewidth = 1, 
              show.legend = FALSE)+
  geom_point(aes(hour,value, color=name),alpha=0.20)+
  facet_wrap(~month)+
  theme_bw()+xlab("Hour")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                'd_can_leaf' = expression("("*T[sim] - T[air]*")")),
                     values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top")+
  NULL

df_res_WR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*100,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("VPD [Pa]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_WR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(Tair,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab("T air [ºC]")+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_WR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric(),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(WS,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("WS [m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL

df_res_WR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(PPFD,d_can_air,color=hour))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("PPFD ["*mu*"mol"*m^{-2}~s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=12)+
  NULL


df_res_WR %>%
  mutate(d_can_air = Tcan-Tair,
         d_can_leaf = tcleaf-Tair,
         hour = lubridate::hour(timestamp),
         month = lubridate::month(timestamp)) %>% 
  ggplot(aes(VPD*WS,d_can_air,color=PPFD))+ 
  geom_abline(slope = 0,intercept=0)+
  geom_point(shape=19)+
  geom_smooth(method="lm", color = "grey30")+
  theme_bw()+xlab(expression("VPD x WS [Pa m "*s^{-1}*"]"))+ylab(expression(Delta*"T [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL


WR_1 <- df_res_WR %>%
  filter(PPFD>0,LE>0) %>% 
  ggplot(aes(lE,LE))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 100, y = 675, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_WR$lE, df_res_WR$LE,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated Latent Heat [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

WR_2 <- df_res_WR %>%
  ggplot(aes(Rnet, NETRAD))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 150, y = 751, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_WR$Rnet, df_res_WR$NETRAD,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  xlab(expression("Simulated NET RADIATION [J "*m^{-2}*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

WR_3 <- df_res_WR %>%
  ggplot(aes(tcleaf, Tcan))+ 
  geom_point(shape=21, color="grey")+
  geom_smooth(method="lm", color = "grey30")+
  geom_abline(slope=1,intercept=0, linetype=2, color="red4",alpha=0.4)+
  annotate("text", x = 5, y = 34, size = 4,
           label = paste("Pearson r = ", 
                         signif(cor(df_res_WR$tcleaf, df_res_WR$Tcan,use="complete.obs"),3)))+
  theme_bw()+ylab(expression("Observed Canopy temperature [ºC]"))+
  xlab(expression("Simulated Leaf temperature [ºC]"))+
  theme(legend.text = element_text(size = 16))+
  scale_color_gradient2(low="cyan3", mid="gold",high="orange",midpoint=1000)+
  NULL

ggpubr::ggarrange(WR_1,WR_2,WR_3, labels = c("a","b","c"),ncol = 3)
#T
plot.eval.dens(df_res_WR$tcleaf,df_res_WR$Tcan)



######################################
<<<<<<< HEAD











=======
library(lme4)
library(lmerTest)
library(relaimpo)
library(DAAG)
library(party)
library(rpart)
library(rpart.plot)
library(mlbench)
library(caret)
library(pROC)
library(tree)
library(ggplot2)
library(ggpubr)

create_train_test <- function(data, size = 0.8) {
  ind <- sample(2, nrow(data), replace = T, prob = c(size, 1-size))
    return (list(data[ind == 1, ],data[ind == 2, ]))
}

accuracy_tune <- function(fit) {
  predict_unseen <- predict(fit, data_test, type = 'class')
  table_mat <- table(data_test$cat_zero, predict_unseen)
  accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
  accuracy_Test
}
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b


#### ALL ####
# set.seed(1234)
# 
# bind_rows(df_res_HA,df_res_ME,df_res_NR,df_res_WR) %>% 
#   mutate(d_can_air = Tcan-Tair,
#          d_can_leaf = tcleaf-Tair,
#          VPD = VPD*100,
#          hour = lubridate::hour(timestamp) %>% as.numeric,
#          month = lubridate::month(timestamp),
#          cat_zero = case_when(d_can_air<0~"<0",
#                               TRUE~'>=0'),
#          cat_zero = as.factor(cat_zero)) %>%
#   dplyr::select(LE,WS,VPD,LAI,PPFD,#FAPAR,
#                 cat_zero,MAP,MAT,z) ->faa1
# 
# data_train1 <- create_train_test(faa1, 0.8, train = TRUE)
# data_test1 <- create_train_test(faa1, 0.8, train = FALSE)
# dim(data_train1)
# 
# fit1 <- rpart(cat_zero~., data = data_train1, method = 'class',maxdepth = 4)
# p <- predict(fit1, data_test1, type = "class")
# rpart.plot(fit1,type = 0,box.palette = "Grays")
# # plotcp(fit)
# # printcp(fit)
# rpart.rules(fit1)
# caret::confusionMatrix(p, reference= data_test1$cat_zero, positive='>=0')
# 



#load files available
#### read the files' paths
filenames.fluxnet<- list.files(path="R/data/final_sites", "*.csv$", full.names=TRUE,recursive = TRUE)

sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")

## REMOVE LAI <0.5, filter Tcan > -50, Tcan <70

#Data from all Canopy temperature sites
list(3,6,7,8,9,10,12,13,19,20,21,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,
  42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
  67,68,69,70,71,72) %>% 
  purrr::map(function(x){
    
    filename <- filenames.fluxnet[x]
    fluxnet_data <- data.table::fread(filename, header = T, 
                                      quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                      colClasses = "numeric", integer64 = "character")
    site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]
    
    data_flx_pre <- fluxnet_data %>% 
      cbind(site = site) %>% 
      mutate(timestamp = lubridate::ymd_hms(DateTime)) %>% 
      left_join(sites_metadata) 
    
<<<<<<< HEAD
    if(unique(data_flx_pre$IGBP) %in% c("ENF","DBF")){
      data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)
      
      return(data_flx_pre)
    }
  }) %>% bind_rows() -> df

# BADM <- list.files("R/data/BADM",full.names = TRUE) %>% 
#   purrr::map(read_excel) %>% bind_rows()
# 
# get_spp_AMF(site = sites$site[3], BADM)
# # Run only first time
# # Define coordinate reference system

# Checking plot
theme_set(theme_bw())
world <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf') 

sites <- df %>% dplyr::select(site,lat,long,IGBP) %>% distinct()

# Create sf object
lat_long_sf <- st_as_sf(sites[,-4], coords = c("long", "lat"), crs = sf::st_crs(world))


world <- st_transform(world, crs = st_crs("ESRI:102003"))
lat_long_sf <- st_transform(lat_long_sf, crs = st_crs("ESRI:102003"))

ggplot() +
  geom_sf(data = world)+
  geom_sf(data=lat_long_sf, aes(color=sites$IGBP),
          # colour="Deep Pink",
          # fill="Pink",
          pch=21, size=3, alpha=I(1))+
  coord_sf(xlim = c(-3500000, 3000000), ylim = c(-1000000, 4000000))+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name=""
  )








gc()
set.seed(1234)

df %>% 
=======
    data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)
    
    return(data_flx_pre)
    
  }) %>% bind_rows() -> df

df %>% dplyr::select(site,IGBP) %>% unique()

df2 <- df %>% 
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
  filter(LAI >= 1, Tcan > 0, Tcan < 70, PPFD>=25, IGBP %in% c("DBF","ENF"#,"MF","EBF"
                                                              )) %>%
  # drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,Tcan, Ustar,LE,SWC) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
<<<<<<< HEAD
                         TRUE~GPP)) %>%
=======
                         TRUE~GPP))
df2 %>% dplyr::select(site,IGBP) %>% unique()

set.seed(1234)

df2 %>%
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         cat_zero = case_when(d_can_air<0~"<0",
                              TRUE~'>=0'),
<<<<<<< HEAD
         cat_zero = as.factor(cat_zero),
         EF=LE/NETRAD) %>%
  dplyr::select(LE,WS,VPD,CO2,IGBP,
                Tair,LAI,PPFD,EF,#MAT,MAP,
                # z,
                Ustar,
                # elev,
                FAPAR,P,rH,SWC,
=======
         cat_zero = as.factor(cat_zero)) %>%
  dplyr::select(LE,WS,VPD,CO2,IGBP,
                Tair,LAI,PPFD,MAT,MAP,z,Ustar,elev,FAPAR,P,rH,SWC,
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
                GPP,#VPD_leaf_calc,
                cat_zero) ->faa2

train_test <- create_train_test(faa2, 0.7)
data_train2 <- train_test[[1]]
data_test2 <- train_test[[2]]
dim(data_train2)

fit2 <- rpart(cat_zero~., data = data_train2, method = 'class',
              control = list(minsplit = 500, xval = 20, maxdepth = 5))
p <- predict(fit2, data_test2, type = "class")
<<<<<<< HEAD
rpart.plot(fit2,type = 2,box.palette = "Grays",branch = 1,cex=1,compress=TRUE)
=======
rpart.plot(fit2,type = 2,box.palette = "Grays")
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
# plotcp(fit2)
# printcp(fit)
# rpart.rules(fit2)
caret::confusionMatrix(p, reference= data_test2$cat_zero, positive='>=0')
<<<<<<< HEAD
# 
# 
# df %>% 
#   filter(LAI >= 1, Tcan > 0, Tcan < 70, PPFD>=25, IGBP %in% c("DBF","ENF"#,"MF","EBF"
#   )) %>%
#   # drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,Tcan, Ustar,LE,SWC) %>% 
#   mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
#          GPP = case_when(GPP<0~0,
#                          TRUE~GPP)) %>%
#   mutate(d_can_air = Tcan-Tair,
#          VPD = VPD*100,
#          hour = lubridate::hour(timestamp) %>% as.numeric,
#          month = lubridate::month(timestamp),
#          cat_zero = case_when(d_can_air<0~"<0",
#                               TRUE~'>=0'),
#          cat_zero = as.factor(cat_zero)) %>%
#   dplyr::select(LE,WS,VPD,CO2,IGBP,
#                 Tair,LAI,PPFD,MAT,MAP,
#                 z,
#                 Ustar,
#                  elev,
#                 FAPAR,P,rH,SWC,
#                 GPP,
#                 d_can_air) ->faa2
# 
# train_test <- create_train_test(faa2, 0.7)
# data_train2 <- train_test[[1]]
# data_test2 <- train_test[[2]]
# dim(data_train2)
# 
# fit2 <- rpart(d_can_air ~., data = data_train2, method = 'anova',
#               control = list(minsplit = 500, xval = 10, maxdepth = 5))
# p <- predict(fit2, data_test2)
# rpart.plot(fit2,type = 2,box.palette = "Grays")
# sqrt(mean((data_test2$d_can_air-p)^2,na.rm=TRUE))
# (cor(data_test2$d_can_air,p, use="complete.obs"))^2













=======
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b







#### Models DeltaT ####
<<<<<<< HEAD
faa3 <- df %>% 
=======
df %>% 
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
  filter(LAI >= 1, Tcan > 0, Tcan < 70, IGBP %in% c("DBF","ENF"#,"MF","EBF"
  )) %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
<<<<<<< HEAD
         VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD),
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         year = lubridate::year(timestamp),
         EF=LE/NETRAD
         ) %>% 
  group_by(site) %>%
  mutate(max_SWC = max(SWC,na.rm=TRUE),
         min_SWC = min(SWC,na.rm=TRUE)) %>%
  rowwise() %>%
  mutate(REW = (max_SWC-SWC)/(max_SWC-min_SWC))

faa3 %>% 
  filter(PPFD>=25) %>% 
  ggplot(aes(EF,d_can_air))+
  geom_point()+
  geom_vline(xintercept=1,linetype=2)+
  geom_vline(xintercept=-1,linetype=2)+
  geom_vline(xintercept=0,linetype=2)+
  xlim(c(-2,2))

# mod_IGBP <- lmer(d_can_air~IGBP+MAP+MAT+elev+z+(1|site),data=faa3 %>% filter(P==0))
# summary(mod_IGBP)
# MuMIn::r.squaredGLMM(mod_IGBP)

mod <- lmer(d_can_air~WS+(PPFD+Tair+rH+EF+LAI)*IGBP+
              (1|site),data=faa3%>% filter(EF>=0,EF<=1))

mod <- lmer(d_can_air~(WS+PPFD+Tair)*IGBP+
              (1|site),data=faa3%>% filter(EF>=0.9,EF<=1))
mod <- lmer(d_can_air~(WS+PPFD+Tair+REW+EF)*IGBP+
              (1|site),data=faa3%>% filter(EF>=0.9,EF<=1))
# mod <- lmer(d_can_air~(LE+WS+NETRAD+Tair+LAI+REW)*IGBP+
#               (1|site),data=faa3)#%>% filter(P==0))
# 
# mod <- lmer(d_can_air ~ LE + WS + NETRAD + Tair + REW + IGBP + (1 | site) + LE:IGBP + WS:IGBP + NETRAD:IGBP + Tair:IGBP + REW:IGBP,data=faa3)
# # 
# mod <- lmer(d_can_air ~ VPD + WS + PPFD + Tair + LE + P + H + IGBP +
#               (1 | site) + (1 | hour) + VPD:IGBP + WS:IGBP + PPFD:IGBP + Tair:IGBP + LE:IGBP + P:IGBP + H:IGBP ,
#             data=faa3)
summary(mod)
step(mod)
MuMIn::r.squaredGLMM(mod)
car::vif(mod)
performance::check_collinearity(mod)


library(effects)
closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)


# 
# eff<-effect("VPD_leaf_calc:IGBP", partial.residuals=T,  mod)
# x.fit <- unlist(eff$x.all$VPD_leaf_calc)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, VPD_leaf_calc = eff$x$VPD_leaf_calc, IGBP=eff$x$IGBP)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$VPD_leaf_calc)] + eff$residuals)
# 
# p1 <- ggplot(x, aes(x = VPD_leaf_calc, y = fit, color=IGBP)) +
#   # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
#   #            shape = 19, show.legend = FALSE, alpha=0.01) +
#   geom_abline(slope=0,intercept=0, linetype=2)+
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), linetype=2) +
#   geom_line(aes(y= upper), linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   theme_bw()+xlab("VPD  leaf [Pa]")+ylab(expression(Delta*"T [ºC]"))+
#   scale_color_manual(values = c("#E69F00", "#56B4E9"),
#                      name =""
#   )+
#   scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#                     name =""
#   )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top",
#         axis.text= element_text(size = 14),
#         axis.title=element_text(size=14))
# 
# 
# 
# eff<-effect("VPD:IGBP", partial.residuals=T,  mod)
# x.fit <- unlist(eff$x.all$VPD)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, VPD = eff$x$VPD, IGBP=eff$x$IGBP)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$VPD)] + eff$residuals)
# 
# p1 <- ggplot(x, aes(x = VPD, y = fit, color=IGBP)) +
#   # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
#   #            shape = 19, show.legend = FALSE, alpha=0.01) +
#   geom_abline(slope=0,intercept=0, linetype=2)+
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), linetype=2) +
#   geom_line(aes(y= upper), linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   theme_bw()+xlab("VPD [Pa]")+ylab(expression(Delta*"T [ºC]"))+
#   scale_color_manual(values = c("#E69F00", "#56B4E9"),
#                      name =""
#   )+
#   scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#                     name =""
#   )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top",
#         axis.text= element_text(size = 14),
#         axis.title=element_text(size=14))




# eff<-effect("WS:IGBP", partial.residuals=T,  mod)
# x.fit <- unlist(eff$x.all$WS)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, WS = eff$x$WS, IGBP=eff$x$IGBP)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$WS)] + eff$residuals)
# 
# p2 <- ggplot(x, aes(x = WS, y = fit, color=IGBP)) +
#   # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
#   #            shape = 19, show.legend = FALSE, alpha=0.01) +
#   geom_abline(slope=0,intercept=0, linetype=2)+
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), linetype=2) +
#   geom_line(aes(y= upper), linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   theme_bw()+xlab("WS [m"~s^{-1}*"]")+ylab(expression(Delta*"T [ºC]"))+
#   scale_color_manual(values = c("#E69F00", "#56B4E9"),
#                      name =""
#   )+
#   scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#                     name =""
#   )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top",
#         axis.text= element_text(size = 14),
#         axis.title=element_text(size=14))


eff<-effect("WS", partial.residuals=T, mod)
x.fit <- unlist(eff$x.all)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, WS = eff$x$WS)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$WS)] + eff$residuals)

p2 <- ggplot(x, aes(x = WS, y = fit)) +
  # geom_point(data = xy, aes(x = x, y = y),
  #            col = "grey60", shape = 1, show.legend = FALSE) +
  geom_abline(slope=0,intercept=0, linetype=2)+
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  theme_bw()+xlab("WS ["*m~s^{-1}*"]")+ylab(expression(Delta*"T [ºC]"))



eff<-effect("PPFD:IGBP", partial.residuals=T,  mod)
x.fit <- unlist(eff$x.all$PPFD)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, PPFD = eff$x$PPFD, IGBP=eff$x$IGBP)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$PPFD)] + eff$residuals)

p3 <- ggplot(x, aes(x = PPFD, y = fit, color=IGBP)) +
  # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
  #            shape = 19, show.legend = FALSE, alpha=0.01) +
  geom_abline(slope=0,intercept=0, linetype=2)+
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  theme_bw()+xlab("PPFD ["*mu*"mol"*m^{-2}~s^{-1}*"]")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))





eff<-effect("Tair:IGBP", partial.residuals=T,  mod)
x.fit <- unlist(eff$x.all$Tair)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, Tair = eff$x$Tair, IGBP=eff$x$IGBP)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$Tair)] + eff$residuals)

p4 <- ggplot(x, aes(x = Tair, y = fit, color=IGBP)) +
  # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
  #            shape = 19, show.legend = FALSE, alpha=0.01) +
  geom_abline(slope=0,intercept=0, linetype=2)+
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  theme_bw()+xlab("Tair [ºC]")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



eff<-effect("rH:IGBP", partial.residuals=T,  mod)
x.fit <- unlist(eff$x.all$rH)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, 
                rH = eff$x$rH, IGBP=eff$x$IGBP)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$rH)] + eff$residuals)

p5 <- ggplot(x, aes(x =rH, y = fit, color=IGBP,fill=IGBP)) +
  # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP),
  #            shape = 19, show.legend = FALSE, alpha=0.01) +
  geom_abline(slope=0,intercept=0, linetype=2)+
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  theme_bw()+xlab("rH [-]")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



eff<-effect("EF:IGBP", partial.residuals=T,  mod)
x.fit <- unlist(eff$x.all$EF)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, EF = eff$x$EF, IGBP=eff$x$IGBP)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$EF)] + eff$residuals)

p1 <- ggplot(x, aes(x = EF, y = fit, color=IGBP,fill=IGBP)) +
  # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP),
  #            shape = 19, show.legend = FALSE, alpha=0.01) +
  geom_abline(slope=0,intercept=0, linetype=2)+
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  theme_bw()+xlab("EF [-]")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



# eff<-effect("REW:IGBP", partial.residuals=T,  mod)
# x.fit <- unlist(eff$x.all$REW)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, REW = eff$x$REW, IGBP=eff$x$IGBP)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$REW)] + eff$residuals)
# 
# p5 <- ggplot(x, aes(x = REW, y = fit, color=IGBP,fill=IGBP)) +
#   # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
#   #            shape = 19, show.legend = FALSE, alpha=0.01) +
#   geom_abline(slope=0,intercept=0, linetype=2)+
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), linetype=2) +
#   geom_line(aes(y= upper), linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   theme_bw()+xlab("REW [-]")+ylab(expression(Delta*"T [ºC]"))+
#   scale_color_manual(values = c("#E69F00", "#56B4E9"),
#                      name =""
#   )+
#   scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#                     name =""
#   )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top",
#         axis.text= element_text(size = 14),
#         axis.title=element_text(size=14))



eff<-effect("LAI:IGBP", partial.residuals=T,  mod)
x.fit <- unlist(eff$x.all$LAI)
trans <- I
x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, LAI = eff$x$LAI, IGBP=eff$x$IGBP)
xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$LAI)] + eff$residuals)

p6 <- ggplot(x, aes(x = LAI, y = fit, color=IGBP,fill=IGBP)) +
  # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP),
  #            shape = 19, show.legend = FALSE, alpha=0.01) +
  geom_abline(slope=0,intercept=0, linetype=2)+
  geom_line(linewidth = 1) +
  geom_line(aes(y= lower), linetype=2) +
  geom_line(aes(y= upper), linetype=2) +
  # geom_smooth(data = xy, aes(x = trans(x), y = y),
  #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
  theme_bw()+xlab("LAI ["*m^2~m^{-2}*"]")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))
# 
# 
# 
# 
# eff<-effect("LAI", partial.residuals=T, mod)
# x.fit <- unlist(eff$x.all)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, LAI = eff$x$LAI)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$LAI)] + eff$residuals)
# 
# p5 <- ggplot(x, aes(x = LAI, y = fit)) +
#   # geom_point(data = xy, aes(x = x, y = y), 
#   #            col = "grey60", shape = 1, show.legend = FALSE) +
#   geom_abline(slope=0,intercept=0, linetype=2)+
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), linetype=2) +
#   geom_line(aes(y= upper), linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   theme_bw()+xlab("LAI ["*m^2~m^{-2}*"]")+ylab(expression(Delta*"T [ºC]"))



# 
# 
# eff<-effect("LE:IGBP", partial.residuals=T,  mod)
# x.fit <- unlist(eff$x.all$LE)
# trans <- I
# x <- data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, LE = eff$x$LE, IGBP=eff$x$IGBP)
# xy <- data.frame(x = x.fit, y = x$fit[closest(trans(x.fit), x$LE)] + eff$residuals)
# 
# p6 <- ggplot(x, aes(x = LE, y = fit, color=IGBP,fill=IGBP)) +
#   # geom_point(data = xy, aes(x = x, y = y, color=eff$x.all$IGBP,fill=eff$x.all$IGBP), 
#   #            shape = 19, show.legend = FALSE, alpha=0.01) +
#   geom_abline(slope=0,intercept=0, linetype=2)+
#   geom_line(linewidth = 1) +
#   geom_line(aes(y= lower), linetype=2) +
#   geom_line(aes(y= upper), linetype=2) +
#   # geom_smooth(data = xy, aes(x = trans(x), y = y),
#   #             method = "loess", span = 2/3, linetype = "dashed", se = FALSE)+
#   theme_bw()+xlab(expression("LE [W "*m^{-2}*"]"))+ylab(expression(Delta*"T [ºC]"))+
#   scale_color_manual(values = c("#E69F00", "#56B4E9"),
#                      name =""
#   )+
#   scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#                     name =""
#   )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top",
#         axis.text= element_text(size = 14),
#         axis.title=element_text(size=14))



mod_plot <- ggarrange(p3,p4,p1,p2,p5,p6,
                      align='hv', labels=c('a', 'b','c',"d","e","f"),
                      ncol=3, nrow = 2, common.legend = T)+
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


ggsave("R/PLOTS/deltaT_model.pdf", plot = mod_plot, width = 36, height = 24, units = "cm")
# ggsave("R/PLOTS/deltaT_model_vpd_leaf.pdf", plot = mod_plot, width = 36, height = 24, units = "cm")




## DeltaT ZERO model ####
coef_DBF <- tibble(intercept = fixef(mod)["(Intercept)"],
                   VPD =  fixef(mod)["VPD"],
                   WS =  fixef(mod)["WS"],
                   PPFD =  fixef(mod)["PPFD"],
                   Tair =  fixef(mod)["Tair"],
                   LAI =  fixef(mod)["LAI"],
                   REW =  fixef(mod)["REW"])

deltaT_zero_DBF <- function(vpd,ws,ppfd,lai,rew){
  foo <- tibble(vpd=vpd,
                ws=ws,
                ppfd=ppfd,
                lai=lai,
                rew=rew,
                Tair=(coef_DBF$intercept+coef_DBF$VPD*vpd+coef_DBF$WS*ws+coef_DBF$PPFD*ppfd+coef_DBF$LAI*lai+coef_DBF$REW*rew)/-coef_DBF$Tair
  )
return(foo)
}

input_deltaT <- expand.grid(seq(0,5000,10),seq(0,2000,1))
p1 <- deltaT_zero_DBF(input_deltaT$Var1,1,input_deltaT$Var2,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>% 
  ggplot(aes(vpd,ppfd))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


input_deltaT <- expand.grid(seq(0,5000,10),seq(0,20,0.1))
deltaT_zero_DBF(input_deltaT$Var1,input_deltaT$Var2,1000,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>% 
  ggplot(aes(vpd,ws))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


input_deltaT <- expand.grid(seq(0,20,0.10),seq(0,2000,1))
deltaT_zero_DBF(1000,input_deltaT$Var1,input_deltaT$Var2,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>% 
  ggplot(aes(ppfd,ws))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



coef_ENF <- tibble(intercept = fixef(mod)["(Intercept)"]+fixef(mod)["IGBPENF"],
                   VPD =  fixef(mod)["VPD"]+fixef(mod)["VPD:IGBPENF"],
                   WS =  fixef(mod)["WS"]+fixef(mod)["WS:IGBPENF"],
                   PPFD =  fixef(mod)["PPFD"]+fixef(mod)["PPFD:IGBPENF"],
                   Tair =  fixef(mod)["Tair"]+fixef(mod)["Tair:IGBPENF"],
                   LAI =  fixef(mod)["LAI"]+fixef(mod)["LAI:IGBPENF"],
                   REW =  fixef(mod)["REW"]+fixef(mod)["REW:IGBPENF"])

deltaT_zero_ENF <- function(vpd,ws,ppfd,lai,rew){
  foo <- tibble(vpd=vpd,
                ws=ws,
                ppfd=ppfd,
                lai=lai,
                rew=rew,
                Tair=(coef_ENF$intercept+coef_ENF$VPD*vpd+coef_ENF$WS*ws+coef_ENF$PPFD*ppfd+coef_ENF$LAI*lai+coef_ENF$REW*rew)/-coef_ENF$Tair
  )
  return(foo)
}


input_deltaT <- expand.grid(seq(0,5000,10),seq(0,2000,1))
p2 <- deltaT_zero_ENF(input_deltaT$Var1,1,input_deltaT$Var2,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>%
  ggplot(aes(vpd,ppfd))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))

input_deltaT <- expand.grid(seq(0,5000,10),10)
deltaT_zero_ENF(input_deltaT$Var1,1,input_deltaT$Var2,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>% 
  ggplot(aes(vpd,Tair))+
  geom_point()+
  geom_point(data=deltaT_negative_ENF(input_deltaT$Var1,10,input_deltaT$Var2,2,0.5),
             aes(vpd,Tair),color="red")+
  # with_blur(
  #   geom_raster(aes(fill = Tair), interpolate = FALSE),
  #   sigma = 5
  # ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


input_deltaT <- expand.grid(seq(0,5000,10),seq(0,20,0.1))
deltaT_zero_ENF(input_deltaT$Var1,input_deltaT$Var2,1000,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>% 
  ggplot(aes(vpd,ws))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))

input_deltaT <- expand.grid(seq(0,20,0.10),seq(0,2000,1))
deltaT_zero_ENF(1000,input_deltaT$Var1,input_deltaT$Var2,2,0.5)%>% 
  filter(Tair>0,Tair<60) %>% 
  ggplot(aes(ppfd,ws))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        # legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))






#### Precipitation effect ####
df %>% 
  filter(LAI >= 1, Tcan > 0, Tcan < 70,PPFD>25, IGBP %in% c("DBF","ENF"#,"MF","EBF"
  )) %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         P_group = case_when(P==0~"[0]",
                             P>0 & P<=20~"(0-20]",
                             TRUE ~ "[>20]"),
         P_group = factor(P_group,
                          labels = c("[0]","(0-20]","[>20]"),
                          levels=c("[0]","(0-20]","[>20]"))) %>% 
ggplot( aes(y = d_can_air, x =P_group, color=IGBP))+
  ggdist::stat_halfeye(
    aes(fill=IGBP),
    color="black",
    adjust = .5, 
    width = .6, 
    .width = c(.5, .95), 
    position =position_dodge(width = 0.7),
    show.legend=FALSE) + 
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .3, 
    ## add some transparency
    alpha = .05,
    position =position_dodge(width = 0.7)
  ) +
  theme_bw()+xlab("Precipitation group [mm"~h^{-1}*"]")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



mod_p <- lmer(d_can_air~P_group*IGBP+(1|site),data=df %>% 
                filter(LAI >= 1, Tcan > 0, Tcan < 70, PPFD>25, IGBP %in% c("DBF","ENF"#,"MF","EBF"
                )) %>%
                mutate(d_can_air = Tcan-Tair,
                       VPD = VPD*100,
                       hour = lubridate::hour(timestamp) %>% as.numeric,
                       month = lubridate::month(timestamp),
                       P_group = case_when(P==0~"[0]",
                                           P>0 & P<=20~"(0-20]",
                                           TRUE ~ "[>20]")))


summary(mod_p)
anova(mod_p)




#### sub daily trends ####
=======
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp))->faa3
mod <- lmer(d_can_air~(VPD+WS+PPFD+elev+LAI+rH+SWC)*IGBP+(1|site),data=faa3)
summary(mod)
MuMIn::r.squaredGLMM(mod)


>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b
#geom_smooth with random effect (https://gist.github.com/richardbeare/b679c38dcb644ec50ea34ac061200504)
predictdf.gam <- function(model, xseq, se, level) {
  olddata <- model.frame(model)
  if (is.null(olddata$randomid)) {
    newdata= tibble(x=xseq)
  } else {
    newdata = tibble(x=xseq, randomid=olddata$randomid[1])
  }
  pred <- predict(model, exclude="s(randomid)", newdata = newdata,
                  se.fit = se, level = level, interval = if (se)
                    "confidence"
                  else "none")
  if (se) {
    y = pred$fit
    ci <- pred$se.fit * 1.96
    ymin = y - ci
    ymax = y + ci
    tibble(x = xseq, y, ymin, ymax, se = pred$se.fit)
  }
  else {
    tibble(x = xseq, y = as.vector(pred))
  }
  
}
environment(predictdf.gam) <- environment(ggplot2:::predictdf.glm)

df %>% 
  filter(LAI >= 1, Tcan > 0, Tcan < 70, IGBP %in% c("DBF","ENF"#,"MF","EBF"
                                                      )) %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp)) %>%
  ggplot() +
  geom_abline(slope=0,intercept=0,linetype=2)+
  # geom_point(aes(hour,d_can_air, color=IGBP),alpha=0.10)+
  # geom_smooth(mapping = aes(hour,d_can_air, color=IGBP, fill=IGBP),linewidth = 1, 
  #             show.legend = TRUE)+
  geom_smooth(aes(hour,d_can_air,
                  group=IGBP, colour=IGBP, fill=IGBP, randomid=factor(site)),
              method="gam", formula = y ~ s(x, bs = "cs") + s(randomid,
                                                              bs="re"))+
  facet_wrap(~month)+
  theme_bw()+xlab("Hour")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top")+
  NULL
<<<<<<< HEAD


p1 <- df %>%
  filter(LAI >= 1, Tcan > 0, Tcan < 70, IGBP %in% c("DBF","ENF"
  )) %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         month_lab = lubridate::month(timestamp,label = T) ,
         day = lubridate::day(timestamp),
         year = lubridate::year(timestamp) %>% as.factor()) %>%
  group_by(month_lab,month,day,IGBP) %>% 
  summarise(max_deltaT = quantile(d_can_air,0.95,na.rm=TRUE),
            min_deltaT = quantile(d_can_air,0.05,na.rm=TRUE)) %>% 
  group_by(month_lab,month,IGBP) %>% 
  summarise(sd_max = sd(max_deltaT,na.rm=TRUE),
            sd_min = sd(min_deltaT,na.rm=TRUE),
            max_deltaT = mean(max_deltaT,na.rm=TRUE),
            min_deltaT = mean(min_deltaT,na.rm=TRUE)) %>% 
  ungroup() %>%
  # mutate(month = lubridate::month(month) %>% as.factor()) %>% 
  ggplot()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_point(aes(month,max_deltaT,color=IGBP),shape=1)+
  geom_errorbar(aes(x=month,ymax=max_deltaT+sd_max,ymin=max_deltaT-sd_max,color=IGBP),
                width=.1)+
  geom_line(aes(month,max_deltaT,color=IGBP))+
  geom_point(aes(month,min_deltaT,color=IGBP),shape=2)+
  geom_errorbar(aes(x=month,ymax=min_deltaT+sd_min,ymin=min_deltaT-sd_min,color=IGBP),
                width=.1)+
  geom_line(aes(month,min_deltaT,color=IGBP), linetype=3)+
  theme_bw()+xlab("Month")+ylab(expression(Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.ticks.length.x = unit(0.5, "cm"),
        axis.text.x = element_text(vjust = 5.5,
                                   hjust = -0.2))+
  scale_x_continuous(breaks = seq(1,12,1), 
                     labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul",
                                "Aug","Sep","Oct","Nov","Dec"))+
  NULL

# 
# 
# df %>%
#   filter(LAI >= 1, Tcan > 0,  Tair > 0, Tcan < 70, IGBP %in% c("DBF","ENF"
#   )) %>%
#   mutate(d_can_air = (Tcan-Tair),
#          VPD = VPD*100,
#          hour = lubridate::hour(timestamp) %>% as.numeric,
#          month = lubridate::month(timestamp),
#          month_lab = lubridate::month(timestamp,label = T) ,
#          day = lubridate::day(timestamp),
#          year = lubridate::year(timestamp) %>% as.factor()) %>%
#   group_by(month_lab,month,day,site,IGBP) %>% 
#   summarise(mean_Tair = mean(Tair,na.rm=TRUE),
#             max_deltaT = quantile(d_can_air,0.95,na.rm=TRUE)/mean_Tair,
#             min_deltaT = quantile(d_can_air,0.05,na.rm=TRUE)/mean_Tair) %>% 
#   group_by(month_lab,month,IGBP) %>% 
#   summarise(sd_max = sd(max_deltaT,na.rm=TRUE),
#             sd_min = sd(min_deltaT,na.rm=TRUE),
#             max_deltaT = mean(max_deltaT,na.rm=TRUE),
#             min_deltaT = mean(min_deltaT,na.rm=TRUE)) %>% 
#   ungroup() %>%
#   # mutate(month = lubridate::month(month) %>% as.factor()) %>% 
#   ggplot()+
#   geom_abline(slope=0,intercept=0,linetype=2)+
#   geom_point(aes(month,max_deltaT,color=IGBP),shape=1)+
#   geom_errorbar(aes(x=month,ymax=max_deltaT+sd_max,ymin=max_deltaT-sd_max,color=IGBP),
#                 width=.1)+
#   geom_line(aes(month,max_deltaT,color=IGBP))+
#   geom_point(aes(month,min_deltaT,color=IGBP),shape=2)+
#   geom_errorbar(aes(x=month,ymax=min_deltaT+sd_max,ymin=min_deltaT-sd_max,color=IGBP),
#                 width=.1)+
#   geom_line(aes(month,min_deltaT,color=IGBP))+
#   theme_bw()+xlab("Month")+ylab(expression(Delta*"T/"*T[air]~"[-]"))+
#   scale_color_manual(values = c("#E69F00", "#56B4E9"),
#                      name =""
#   )+
#   scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#                     name =""
#   )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top",
#         axis.ticks.length.x = unit(0.5, "cm"),
#         axis.text.x = element_text(vjust = 5.5,
#                                    hjust = -0.2))+
#   scale_x_continuous(breaks = seq(1,12,1), 
#                      labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul",
#                                 "Aug","Sep","Oct","Nov","Dec"))+
#   NULL


p2 <- df %>%
  filter(LAI >= 1, Tcan > 0, Tcan < 70, IGBP %in% c("DBF","ENF"
  )) %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         month_lab = lubridate::month(timestamp,label = T) ,
         day = lubridate::day(timestamp),
         year = lubridate::year(timestamp) %>% as.factor()) %>%
  group_by(month_lab,month,day,IGBP) %>% 
  summarise(max_deltaT = quantile(d_can_air,0.95,na.rm=TRUE),
            min_deltaT = quantile(d_can_air,0.05,na.rm=TRUE)) %>% 
  group_by(month_lab,month,IGBP) %>% 
  summarise(sd_max = sd(max_deltaT,na.rm=TRUE),
            sd_min = sd(min_deltaT,na.rm=TRUE),
            max_deltaT = mean(max_deltaT,na.rm=TRUE),
            min_deltaT = mean(min_deltaT,na.rm=TRUE)) %>% 
  ungroup() %>%
  ggplot()+
  # geom_abline(slope=0,intercept=0,linetype=2)+
  geom_point(aes(month,sd_max,color=IGBP),shape=1)+
  geom_line(aes(month,sd_max,color=IGBP))+
  geom_point(aes(month,sd_min,color=IGBP),shape=2)+
  geom_line(aes(month,sd_min,color=IGBP),linetype=3)+
  theme_bw()+xlab("Month")+ylab(expression("SD"~Delta*"T [ºC]"))+
  scale_color_manual(values = c("#E69F00", "#56B4E9"),
                     name =""
  )+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    name =""
  )+
  theme(legend.text = element_text(size = 16),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        axis.ticks.length.x = unit(0.5, "cm"),
        axis.text.x = element_text(vjust = 5.5,
                                   hjust = -0.2))+
  scale_x_continuous(breaks = seq(1,12,1), 
                     labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul",
                                "Aug","Sep","Oct","Nov","Dec"))+
  NULL


ggarrange(p1,p2,
          align='hv', labels=c('a', 'b'),
          ncol=1, nrow = 2, common.legend = T)+
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))

# 
# df %>% 
#   filter(LAI >= 1, Tcan > 0, Tcan < 70, IGBP %in% c("DBF"
#   ), !site%in%c('US-CMW', 'US-Ha1','US-Rpf')
#   ) %>%
#   mutate(d_can_air = Tcan-Tair,
#          VPD = VPD*100,
#          hour = lubridate::hour(timestamp) %>% as.numeric,
#          month = lubridate::month(timestamp),
#          year = lubridate::year(timestamp) %>% as.factor()) %>%
#   ggplot() +
#   geom_abline(slope=0,intercept=0,linetype=2)+
#   # geom_point(aes(hour,d_can_air, color=year),alpha=0.10)+
#   # geom_smooth(mapping = aes(hour,d_can_air, color=IGBP, fill=IGBP),linewidth = 1, 
#   #             show.legend = TRUE)+
#   geom_smooth(aes(hour,d_can_air,
#                   group=year, colour=year, fill=year, randomid=factor(site)),
#               method="gam", formula = y ~ s(x, bs = "cs") + s(randomid,bs="re")
#               )+
#   facet_wrap(~month)+
#   theme_bw()+xlab("Hour")+ylab(expression(Delta*"T [ºC]"))+
#   # scale_color_manual(values = c("#E69F00", "#56B4E9"),
#   #                    name =""
#   # )+
#   # scale_fill_manual(values = c("#E69F00", "#56B4E9"),
#   #                   name =""
#   # )+
#   theme(legend.text = element_text(size = 16),
#         strip.background = element_rect(fill="white"),
#         legend.position="top")+
#   NULL


################################################################33

x <- seq(273.15,330,1)
y <- exp(34.494-4924.99/((x-273.15)+237.1))/(((x-273.15)+105)^1.57)

es_mod <- lm(y ~ x + I(x^2) + I(x^3)) 
summary(es_mod)

es = es_mod$coefficients[1]+
  es_mod$coefficients[2]*x+
  es_mod$coefficients[3]*x^2+
  es_mod$coefficients[4]*x^3

plot(x,y)
lines(x,es,col="red")



calc_Tair_DeltaT <- function(Rs = 800, Rl = 400, gs = 0.005, 
                             rH = 100, alpha = 0.5, emi = 0.97,
                             h = 0.5, n=1,
                             v1 = es_mod$coefficients[4] %>% as.numeric(),
                             v2 = es_mod$coefficients[3]%>% as.numeric(),
                             v3 = es_mod$coefficients[2]%>% as.numeric(),
                             v4 = es_mod$coefficients[1]%>% as.numeric()){
  
  sigma = 5.67e-8 #W m-2 K-4
  rho = 1225 #g m-3
  Cp = 1.01 #J g-1 K-1
  gamma = 67 #Pa K-1
  
  
  
  a = emi*sigma/(rho*Cp)
  b = v1 * n*gs*(1-h)/(gamma+gamma*gs*rH)
  c = v2 * n*gs*(1-h)/(gamma+gamma*gs*rH)
  d = v3 * n*gs*(1-h)/(gamma+gamma*gs*rH)
  f = (v4 * n*gs*(1-h)/(gamma+gamma*gs*rH)) - ((emi*Rl+alpha*Rs)/(2*rho*Cp))
  
  res <- tryCatch({
    
    uniroot(function(x){
      # a*x^4 + b*x^3  + c*x^2  + d*x +f
      es <- exp(34.494-4924.99/((x-273.15)+237.1))/(((x-273.15)+105)^1.57)
      ((emi*Rl+alpha*Rs-(2*emi*sigma*x^4))/(2*rho*Cp))-(n*gs*es*(1-h)/(gamma+gamma*gs*rH))
    },lower=273.15,upper=400)->foo
    
    return(foo$root)
  },
  error=function(e) {
    # Choose a return value in case of error
    return(NA)
  }
  )
  
  return(res)
  
}


# calc_Tair_DeltaT()
expand.grid(seq(500,1200,1),seq(0,0.0007,0.00001)) %>% t() %>% as.data.frame() %>%  as.list() %>% 
  purrr::map(function(x){
    data.frame(Tair = calc_Tair_DeltaT(Rs=x[1],Rl=x[1]/2,rH = 100,gs=x[2],h=0.4),
               Rs = x[1],
               gs = x[2])}) %>% 
  bind_rows() -> res_T
# T_lim = (res_T$Rs-500)/10+273.15
# mod_ = lm(gs~Rs*Tair, data=res_T )
# pred_gs_T = predict(mod_, data.frame(Rs = res_T$Rs,Tair=T_lim) %>% distinct()) 
# db = data.frame(Rs = res_T$Rs %>% unique(), gs = pred_gs_T) %>% dplyr::distinct()
T1 <- res_T %>% 
  filter(Tair>=273.15, Tair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  # geom_line(data=db,aes(Rs,gs),se=FALSE)+
  annotate("text", x = 1050, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 100 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 1050, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 40 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))




expand.grid(seq(500,1200,1),seq(0,0.0007,0.00001)) %>% t() %>% as.data.frame() %>%  as.list() %>% 
  purrr::map(function(x){
    data.frame(Tair = calc_Tair_DeltaT(Rs=x[1],Rl=x[1]/2,rH = 1,gs=x[2],h=0.4),
               Rs = x[1],
               gs = x[2])}) %>% 
  bind_rows() -> res_T2
T2 <- res_T2 %>% 
  filter(Tair>=273.15, Tair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 1050, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 1 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 1050, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 40 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



expand.grid(seq(500,1200,1),seq(0,0.0007,0.00001)) %>% t() %>% as.data.frame() %>%  as.list() %>% 
  purrr::map(function(x){
    data.frame(Tair = calc_Tair_DeltaT(Rs=x[1],Rl=x[1]/2,rH = 100,gs=x[2],h=0.9),
               Rs = x[1],
               gs = x[2])}) %>% 
  bind_rows() -> res_T3
T3 <- res_T3 %>% 
  filter(Tair>=273.15, Tair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 1050, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 100 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 1050, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 90 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))



expand.grid(seq(500,1200,1),seq(0,0.0007,0.00001)) %>% t() %>% as.data.frame() %>%  as.list() %>% 
  purrr::map(function(x){
    data.frame(Tair = calc_Tair_DeltaT(Rs=x[1],Rl=x[1]/2,rH = 1,gs=x[2],h=0.9),
               Rs = x[1],
               gs = x[2])}) %>% 
  bind_rows() -> res_T4
T4 <- res_T4 %>% 
  filter(Tair>=273.15, Tair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 1050, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 1 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 1050, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 90 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


ggarrange(T1,T2,T3,T4,labels = c("a","b","c","d"),ncol = 2,nrow=2,
          common.legend = TRUE,legend ="top" )








#### Models DeltaT ####
faa3 <- df %>% 
  filter(LAI >= 1, Tcan > 0, Tcan < 70, LE >=0,
         IGBP %in% c("DBF","ENF"#,"MF","EBF"
  )) %>%
  mutate(d_can_air = Tcan-Tair,
         Tkair = Tair + 273.15,
         VPD = VPD*100,
         VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD),
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         year = lubridate::year(timestamp),
         EF=LE/NETRAD,
         Rs = SW_IN,
         Rs = case_when(is.na(Rs) & !is.na(PPFD) ~ bigleaf::PPFD.to.Rg(PPFD),
                        TRUE ~ Rs),
         leafwidth = case_when(IGBP == 'DBF'~0.05,
                               TRUE ~ 0.002),
         gb = case_when(IGBP == 'DBF'~ 0.007*sqrt(WS/leafwidth),
                        TRUE ~ 0.004*(WS^0.6/leafwidth^0.4)),
         Rb = 1/gb,
         gh = 0.92*gb,
         Rh = 1/gh,
         E = LE.to.ET(LE,Tair)/1000,
         gbs = (PA*1000*E)/VPD_leaf_calc,
         r_s = (1/gbs)-Rb,
         gs=1/r_s
 ) %>% 
  filter(gs>=0,!is.na(Rs),!is.infinite(Rh))
  # group_by(site) %>%
  # mutate(max_SWC = max(SWC,na.rm=TRUE),
  #        min_SWC = min(SWC,na.rm=TRUE)) %>%
  # rowwise() %>%
  # mutate(#REW = (max_SWC-SWC)/(max_SWC-min_SWC),
  #        ,
  #        )

# faa3 %>% 
#   filter(PPFD>=25) %>% 
#   ggplot(aes(EF,d_can_air))+
#   geom_point()+
#   geom_vline(xintercept=1,linetype=2)+
#   geom_vline(xintercept=-1,linetype=2)+
#   geom_vline(xintercept=0,linetype=2)+
#   xlim(c(-2,2))

# mod_IGBP <- lmer(d_can_air~IGBP+MAP+MAT+elev+z+(1|site),data=faa3 %>% filter(P==0))
# summary(mod_IGBP)
# MuMIn::r.squaredGLMM(mod_IGBP)

# mod <- lmer(d_can_air~Rh+(gs + Rs + Tkair + rH )*IGBP+
#               (1|site),data=faa3)
# 
# summary(mod)
# # step(mod)
# MuMIn::r.squaredGLMM(mod)
# car::vif(mod)
# performance::check_collinearity(mod)
# 
# 
# ## DeltaT ZERO model ####
# coef_DBF <- tibble(intercept = fixef(mod)["(Intercept)"],
#                    gs =  fixef(mod)["gs"],
#                    Rs =  fixef(mod)["Rs"],
#                    Rh =  fixef(mod)["Rh"],
#                    rH =  fixef(mod)["rH"],
#                    Tkair =  fixef(mod)["Tkair"])
# 
# deltaT_zero_DBF <- function(gs,Rs,Rh,rH,Tkair){
#   foo <- tibble(gs=gs,
#                 Rs=Rs,
#                 Rh=Rh,
#                 rH=rH,
#                 Tkair=(coef_DBF$intercept+coef_DBF$gs*gs+coef_DBF$Rs*Rs+
#                         coef_DBF$Rh*Rh+coef_DBF$rH*rH)/-coef_DBF$Tkair
#   )
#   return(foo)
# }

set.seed(1234)
# train_test <- create_train_test(faa3 %>% 
#                                   dplyr::select(Tkair,gs,d_can_air,Rs,rH,Rh,MAT,MAP) %>% 
#                                   na.omit(), 
#                                 0.7)
# data_train <- train_test[[1]]
# data_test <- train_test[[2]]
# 
# 
# # -----------------------------------------------------------------------------
# # STEP 1: SELECT TUNING PARAMETERS
# 
# # part a: set range of tuning parameters (layer size and weight decay)
# tune_grid_neural <- expand.grid(size = c(1:5, 10),
#                                 decay = c(0, 0.05, 0.1, 1, 2))
# 
# 
# # part b: set some other consrains to be imposed on network (to keep computation manageable)
# # see: p. 361 of Kuhn & Johnson (2013, 
# max_size_neaural <- max(tune_grid_neural$size)
# max_weights_neural <- max_size_neaural*(nrow(data_train) + 1) + max_size_neaural + 1
# 
# # -----------------------------------------------------------------------------
# # STEP 2: SELECT TUNING METHOD
# # set up train control object, which specifies training/testing technique
# train_control_neural <- trainControl(method = "repeatedcv",
#                                      number = 5,
#                                      p = 0.70)
# 
# # start timer
# start_time <- Sys.time()
# 
# # -----------------------------------------------------------------------------
# # STEP 3: TRAIN MODEL
# 
# # use caret "train" function to train svm
# model_ex <- 
#   train(form = Tkair ~ gs+d_can_air+Rs+rH+Rh+MAT+MAP,
#         data = data_train,
#         method = "nnet",
#         tuneGrid = tune_grid_neural,
#         trControl = train_control_neural,
#         metric = "RMSE", # how to select among models
#         trace = FALSE,
#         maxit = 100,
#         MaxNWts = max_weights_neural
#         ) # don't print output along the way
# 
# # end timer
# total_time <- Sys.time() - start_time
# 
# model_ex$value
# 
# pred_mod <- predict(model_ex,data_test)
# 
# plot(pred_mod,data_test$Tkair)


# model_nn <- sfn_predict_nn(faa3 %>% 
#                  dplyr::select(Tkair,gs,d_can_air,Rs,rH,Rh,MAT,MAP) %>% 
#                  na.omit(),c('gs','d_can_air','Rs','rH','Rh','MAT','MAP'),
#                "Tkair"
#               )

# saveRDS(model_nn,file="R/model_nn.rds")
model_nn <- readRDS(file="R/model_nn.rds")
input_deltaT <- expand.grid(seq(0,1200,1),seq(0,0.0007,0.000001))

Tkair <- predict(model_nn$nn, 
        newdata=data.frame(gs=input_deltaT$Var2,d_can_air=0,
                           Rs=input_deltaT$Var1,rH=40,Rh=100,MAT=20,MAP=500))
p1 <- data.frame(gs=input_deltaT$Var2,Rs=input_deltaT$Var1,Tkair = Tkair) %>% 
  filter(Tkair>=273.15, Tkair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tkair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 900, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 100 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 900, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 40 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


Tkair <- predict(model_nn$nn, 
                 newdata=data.frame(gs=input_deltaT$Var2,d_can_air=0,
                                    Rs=input_deltaT$Var1,rH=40,Rh=1,MAT=20,MAP=500))
p2 <-  data.frame(gs=input_deltaT$Var2,Rs=input_deltaT$Var1,Tkair = Tkair) %>% 
  filter(Tkair>=273.15, Tkair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tkair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 900, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 1 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 900, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 40 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


Tkair <- predict(model_nn$nn, 
                 newdata=data.frame(gs=input_deltaT$Var2,d_can_air=0,
                                    Rs=input_deltaT$Var1,rH=90,Rh=100,MAT=20,MAP=500))
p3 <- data.frame(gs=input_deltaT$Var2,Rs=input_deltaT$Var1,Tkair = Tkair) %>% 
  filter(Tkair>=273.15, Tkair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tkair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 900, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 100 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 900, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 90 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


Tkair <- predict(model_nn$nn, 
                 newdata=data.frame(gs=input_deltaT$Var2,d_can_air=0,
                                    Rs=input_deltaT$Var1,rH=90,Rh=1,MAT=20,MAP=500))
p4 <- data.frame(gs=input_deltaT$Var2,Rs=input_deltaT$Var1,Tkair = Tkair) %>% 
  filter(Tkair>=273.15, Tkair<(273.15+60)) %>% 
  ggplot(aes(Rs,gs))+
  with_blur(
    geom_raster(aes(fill = Tkair), interpolate = FALSE),
    sigma = 5
  ) +
  annotate("text", x = 900, y = 0.00015,hjust = 0,
           label = "paste(r[h],\" = 1 s \", m ^ {-1})",parse=TRUE)+
  annotate("text", x = 900, y = 0.0001,hjust = 0,
           label = "paste(H,\" = 90 %\")",parse=TRUE)+
  scale_fill_viridis_c(option = "inferno")+
  theme_bw()+xlab(expression(R[s]~"[W "*m^{-2}*"]"))+ylab(expression(g[s]~"[m "*s^{-1}*"]"))+
  theme(legend.text = element_text(size = 12),
        strip.background = element_rect(fill="white"),
        legend.position="top",
        legend.key.width=unit(5, "cm"),
        axis.text= element_text(size = 14),
        axis.title=element_text(size=14))


ggarrange(p1,p2,p3,p4,labels = c("a","b","c","d"),ncol = 2,nrow=2,
          common.legend = TRUE,legend ="top" )
=======
>>>>>>> aada352b0b74d0b4832b24ee9e2235367baf534b

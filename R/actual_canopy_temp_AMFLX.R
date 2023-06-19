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
library(ggplot2)
library(viridis)
library(Metrics)
library(sf)
library(R.utils)
library(rsplashtest)
source("R/evaluation_sim_plot.R")
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




precip_ERA5 <- read_csv(path="R/data/ERA_values_prec_amf.csv$")%>% 
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
    
    data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)
    
    return(data_flx_pre)
    
  }) %>% bind_rows() -> df

df %>% dplyr::select(site,IGBP) %>% unique()

df2 <- df %>% 
  filter(LAI >= 1, Tcan > 0, Tcan < 70, PPFD>=25, IGBP %in% c("DBF","ENF"#,"MF","EBF"
                                                              )) %>%
  # drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,Tcan, Ustar,LE,SWC) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP))
df2 %>% dplyr::select(site,IGBP) %>% unique()

set.seed(1234)

df2 %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp),
         cat_zero = case_when(d_can_air<0~"<0",
                              TRUE~'>=0'),
         cat_zero = as.factor(cat_zero)) %>%
  dplyr::select(LE,WS,VPD,CO2,IGBP,
                Tair,LAI,PPFD,MAT,MAP,z,Ustar,elev,FAPAR,P,rH,SWC,
                GPP,#VPD_leaf_calc,
                cat_zero) ->faa2

train_test <- create_train_test(faa2, 0.7)
data_train2 <- train_test[[1]]
data_test2 <- train_test[[2]]
dim(data_train2)

fit2 <- rpart(cat_zero~., data = data_train2, method = 'class',
              control = list(minsplit = 500, xval = 20, maxdepth = 5))
p <- predict(fit2, data_test2, type = "class")
rpart.plot(fit2,type = 2,box.palette = "Grays")
# plotcp(fit2)
# printcp(fit)
# rpart.rules(fit2)
caret::confusionMatrix(p, reference= data_test2$cat_zero, positive='>=0')







#### Models DeltaT ####
df %>% 
  filter(LAI >= 1, Tcan > 0, Tcan < 70, IGBP %in% c("DBF","ENF"#,"MF","EBF"
  )) %>%
  mutate(d_can_air = Tcan-Tair,
         VPD = VPD*100,
         hour = lubridate::hour(timestamp) %>% as.numeric,
         month = lubridate::month(timestamp))->faa3
mod <- lmer(d_can_air~(VPD+WS+PPFD+elev+LAI+rH+SWC)*IGBP+(1|site),data=faa3)
summary(mod)
MuMIn::r.squaredGLMM(mod)


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

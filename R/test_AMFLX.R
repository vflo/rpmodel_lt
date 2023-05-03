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
source("R/evaluation_sim_plot.R")
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
########################################################################
# # Reading Eddy Covariance data using the readFlux Package
########################################################################
########################################################################
# 01.find the filename
########################################################################
#load files available
#### read the files' paths
filenames.fluxnet<- list.files(path="R/data/final_sites", "*.csv$", full.names=TRUE,recursive = TRUE)

########################################################################
# 02.Read the data US-Ha
########################################################################
filename <- filenames.fluxnet[10]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

#### read metadata ####
sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

data_flx <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP))
summary(data_flx)

res_HA <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                           co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                           ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA,# data_flx$Ustar, 
                           canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                           elv = data_flx$elev, z = data_flx$z, leafwidth = 0.05, 
                           netrad = NA, #data_flx$NETRAD,
                           beta = 146.0, c_cost = 0.41, 
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, 
                           upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE,
                           energy_params = list(
                             epsleaf = 0.97, #thermal absorptivity of the leaf
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_HA <- as_tibble(res_HA) %>% cbind(data_flx)

res2_HA <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                            co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                            ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                            canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                            elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001,
                            netrad = NA,#data_flx$NETRAD, 
                            beta = 146.0, c_cost = 0.41,
                            do_leaftemp = FALSE,  gb_method = "Choudhury_1988",
                            do_acclimation = TRUE, 
                            upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                            acclim_days = 15, weighted_accl = TRUE,
                            energy_params = list(
                              epsleaf = 0.98, #thermal absorptivity of the leaf
                              ste_bolz = 5.67e-8, #W m^-2 K^-4
                              cpm = 75.38, #J mol^-1 ºC-1
                              J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                              frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                              fanir = 0.35 #Fraction of NIR absorbed
                            ))

df_res2_HA <- as_tibble(res2_HA) %>% cbind(data_flx)



range_plot <- c(2510:2581)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_HA$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_HA$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####
# (ha_no <- ggplot() +
#   geom_line(data = df_res2_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,GPP), color = "black")+
#   geom_line(data = df_res2_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#   theme_bw()+xlab("Date")+ylab("GPP")+
#   NULL)
# 
# (ha_yes <- ggplot() +
#   geom_line(data = df_res2_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,GPP), color = "black")+
#   geom_line(data = df_res_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#   theme_bw()+xlab("Date")+ylab("GPP")+
#   NULL)

# (ha_t <- ggplot() +
#   geom_line(data = df_res_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
#   geom_line(data = df_res_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
#   geom_line(data = df_res_HA %>%
#               slice(range_plot),
#             mapping = aes(timestamp,Tcan), color = "green2",linewidth = 1 )+
#   geom_line(data = df_res_HA %>%
#               slice(range_plot) %>% 
#               mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
#             mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
#   geom_line(data = df_res_HA%>%
#               slice(range_plot) %>% 
#               mutate(T_H = H/(gb*0.92*75.38)+Tair),
#             mapping = aes(timestamp,T_H), color = "grey40",linewidth = 1,linetype = 1 )+
#   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#   theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
#   NULL)

(ha_t <-  df_res_HA %>%
   slice(range_plot) %>% 
   mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15,
          T_H = H/(gb*0.92*75.38)+Tair,
          d_can_air = Tcan-Tair,
          d_can_leaf = Tcan-tcleaf,
          d_lw_leaf = T_lw-tcleaf,
          d_h_leaf = T_H-tcleaf) %>% 
   pivot_longer(c(d_can_air,d_can_leaf,d_lw_leaf,d_h_leaf)) %>%
   ggplot() +
   geom_abline(slope=0,intercept=0,linetype=2)+
   geom_line(mapping = aes(timestamp,value, color=name),linewidth = 1 )+
   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
   theme_bw()+xlab("Date")+ylab(expression(Delta*"T [ºC]"))+
   scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                 'd_can_leaf' = expression("("*T[can] - T[leaf]*")"),
                                 'd_lw_leaf' = expression("("*T[LW] - T[leaf]*")"),
                                 'd_h_leaf' = expression("("*T[H] - T[leaf]*")")),
                      values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                      name =""
   )+
   theme(legend.text = element_text(size = 16))+
   NULL)

# df_res_HA_filter <- df_res_HA %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res_HA_filter$assim,df_res_HA_filter$GPP)
# df_res2_HA_filter <- df_res2_HA %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res2_HA_filter$assim,df_res2_HA_filter$GPP)

#H
# plot.eval.dens(df_res_HA$Qc,df_res_HA$H)
plot.eval.dens(df_res_HA$tcleaf,
               (df_res_HA$H/(df_res_HA$gb*0.92*75.38)+df_res_HA$Tair))

#T LW
plot.eval.dens(df_res_HA$tcleaf,(df_res_HA$LW_OUT/(0.96*5.67e-8))^0.25-273.15)


#NETRAD
# plot.eval.dens(df_res_HA$Rnet,df_res_HA$NETRAD)

#T
plot.eval.dens(df_res_HA$tcleaf,df_res_HA$Tcan)

#T
# plot.eval.dens(df_res_HA$Tair,df_res_HA$Tcan)

#USTAR
# plot.eval.dens(df_res_HA$ust,df_res_HA$Ustar)

#LE
# df_res_HA_filter <- df_res_HA %>% filter(SW_IN>25)
# plot.eval.dens(df_res_HA_filter$e*18.01528*2230*1e-6,df_res_HA_filter$LE)
# df_res2_HA_filter <- df_res2_HA %>% filter(SW_IN>25)
# plot.eval.dens(df_res2_HA_filter$e*18.01528*2230*1e-6,df_res2_HA_filter$LE)








########################################################################
# 02.Read the data US-Me2
########################################################################
filename <- filenames.fluxnet[12]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

#### read metadata ####
sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

data_flx <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP))
summary(data_flx)

res_ME <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                           co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                           ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA,# data_flx$Ustar, 
                           canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                           elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001, 
                           netrad = NA, #data_flx$NETRAD,
                           beta = 146.0, c_cost = 0.41, 
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, 
                           upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE,
                           energy_params = list(
                             epsleaf = 0.96, #thermal absorptivity of the leaf
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_ME <- as_tibble(res_ME) %>% cbind(data_flx)

res2_ME <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                            co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                            ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                            canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                            elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001,
                            netrad = NA,#data_flx$NETRAD, 
                            beta = 146.0, c_cost = 0.41,
                            do_leaftemp = FALSE,  gb_method = "Su_2001",
                            do_acclimation = TRUE, 
                            upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                            acclim_days = 15, weighted_accl = TRUE,
                            energy_params = list(
                              epsleaf = 0.98, #thermal absorptivity of the leaf
                              ste_bolz = 5.67e-8, #W m^-2 K^-4
                              cpm = 75.38, #J mol^-1 ºC-1
                              J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                              frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                              fanir = 0.35 #Fraction of NIR absorbed
                            ))

df_res2_ME <- as_tibble(res2_ME) %>% cbind(data_flx)



range_plot <- c(5472:5626)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_ME$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_ME$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

# (me_no <- ggplot() +
#   geom_line(data = df_res2_ME %>%
#               slice(range_plot),
#             mapping = aes(timestamp,GPP), color = "black")+
#   geom_line(data = df_res2_ME %>%
#               slice(range_plot),
#             mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#   theme_bw()+xlab("Date")+ylab("GPP")+
#   NULL)
# 
# (me_yes <- ggplot() +
#   
#   geom_line(data = df_res2_ME %>%
#               slice(range_plot),
#             mapping = aes(timestamp,GPP), color = "black")+
#   geom_line(data = df_res_ME %>%
#               slice(range_plot),
#             mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#   theme_bw()+xlab("Date")+ylab("GPP")+
#   NULL)
# 
# (me_t <- ggplot() +
#     geom_line(data = df_res_ME %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
#     geom_line(data = df_res_ME %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
#     geom_line(data = df_res_ME %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,Tcan), color = "green2",linewidth = 1 )+
#     geom_line(data = df_res_ME %>%
#                 slice(range_plot) %>% 
#                 mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
#               mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
#     geom_line(data = df_res_ME%>%
#                 slice(range_plot) %>% 
#                 mutate(T_H = H/(gb*0.92*75.38)+Tair),
#               mapping = aes(timestamp,T_H), color = "grey40",linewidth = 1,linetype = 1 )+
#     geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#     theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
#     NULL)

(me_t <-  df_res_ME %>%
   slice(range_plot) %>% 
   mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15,
          T_H = H/(gb*0.92*75.38)+Tair,
          d_can_air = Tcan-Tair,
          d_can_leaf = Tcan-tcleaf,
          d_lw_leaf = T_lw-tcleaf,
          d_h_leaf = T_H-tcleaf) %>% 
   pivot_longer(c(d_can_air,d_can_leaf,d_lw_leaf,d_h_leaf)) %>%
   ggplot() +
   geom_abline(slope=0,intercept=0,linetype=2)+
   geom_line(mapping = aes(timestamp,value, color=name),linewidth = 1 )+
   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
   theme_bw()+xlab("Date")+ylab(expression(Delta*"T [ºC]"))+
   scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                 'd_can_leaf' = expression("("*T[can] - T[leaf]*")"),
                                 'd_lw_leaf' = expression("("*T[LW] - T[leaf]*")"),
                                 'd_h_leaf' = expression("("*T[H] - T[leaf]*")")),
                      values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                      name =""
   )+
   theme(legend.text = element_text(size = 16))+
   NULL)

# df_res_ME_filter <- df_res_ME %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res_ME_filter$assim,df_res_ME_filter$GPP)
# df_res2_ME_filter <- df_res2_ME %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res2_ME_filter$assim,df_res2_ME_filter$GPP)

#H
# plot.eval.dens(df_res_ME$Qc,df_res_ME$H)
df_res_ME_filter <- df_res_ME %>% mutate(gb = case_when(gb==0 ~ 1e-1,
                                                               TRUE~gb))
plot.eval.dens(df_res_ME_filter$tcleaf,
               (df_res_ME_filter$H/(df_res_ME_filter$gb*0.92*75.38)+df_res_ME_filter$Tair))

#T LW
plot.eval.dens(df_res_ME$tcleaf,(df_res_ME$LW_OUT/(0.96*5.67e-8))^0.25-273.15)

#NETRAD
# plot.eval.dens(df_res_ME$Rnet,df_res_ME$NETRAD)

#T
plot.eval.dens(df_res_ME$tcleaf,df_res_ME$Tcan)

#T
# plot.eval.dens(df_res_ME$Tair,df_res_ME$Tcan)

#USTAR
# plot.eval.dens(df_res_ME$ust,df_res_ME$Ustar)

#LE
# df_res_ME_filter <- df_res_ME %>% filter(SW_IN>25)
# plot.eval.dens(df_res_ME_filter$e*18.01528*2230*1e-6,df_res_ME_filter$LE)
# df_res2_ME_filter <- df_res2_ME %>% filter(SW_IN>25)
# plot.eval.dens(df_res2_ME_filter$e*18.01528*2230*1e-6,df_res2_ME_filter$LE)










########################################################################
# 02.Read the data US-NR1
########################################################################
filename <- filenames.fluxnet[13]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

#### read metadata ####
sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

data_flx <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP))
summary(data_flx)

res_NR <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                           co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                           ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA,# data_flx$Ustar, 
                           canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                           elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001, 
                           netrad = NA, #data_flx$NETRAD,
                           beta = 146.0, c_cost = 0.41, 
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, 
                           upscaling_method = "daily", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE,
                           energy_params = list(
                             epsleaf = 0.96, #thermal absorptivity of the leaf
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_NR <- as_tibble(res_NR) %>% cbind(data_flx)

res2_NR <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                            co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                            ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                            canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                            elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001,
                            netrad = NA,#data_flx$NETRAD, 
                            beta = 146.0, c_cost = 0.41,
                            do_leaftemp = FALSE,  gb_method = "Su_2001",
                            do_acclimation = TRUE, 
                            upscaling_method = "daily", hour_reference_T = c(10,12,14), 
                            acclim_days = 15, weighted_accl = TRUE,
                            energy_params = list(
                              epsleaf = 0.98, #thermal absorptivity of the leaf
                              ste_bolz = 5.67e-8, #W m^-2 K^-4
                              cpm = 75.38, #J mol^-1 ºC-1
                              J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                              frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                              fanir = 0.35 #Fraction of NIR absorbed
                            ))

df_res2_NR <- as_tibble(res2_NR) %>% cbind(data_flx)



range_plot <- c(503:654)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_NR$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_NR$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

# (nr_no <- ggplot() +
#    geom_line(data = df_res2_NR %>%
#                slice(range_plot),
#              mapping = aes(timestamp,GPP), color = "black")+
#    geom_line(data = df_res2_NR %>%
#                slice(range_plot),
#              mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#    geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#    theme_bw()+xlab("Date")+ylab("GPP")+
#    NULL)
# 
# (nr_yes <- ggplot() +
#     
#     geom_line(data = df_res2_NR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,GPP), color = "black")+
#     geom_line(data = df_res_NR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#     geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#     theme_bw()+xlab("Date")+ylab("GPP")+
#     NULL)
# 
# (nr_t <- ggplot() +
#     geom_line(data = df_res_NR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
#     geom_line(data = df_res_NR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
#     geom_line(data = df_res_NR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,Tcan), color = "green2",linewidth = 1 )+
#     geom_line(data = df_res_NR %>%
#                 slice(range_plot) %>% 
#                 mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
#               mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
#     geom_line(data = df_res_NR%>%
#                 slice(range_plot) %>% 
#                 mutate(T_H = H/(gb*0.92*75.38)+Tair),
#               mapping = aes(timestamp,T_H), color = "grey40",linewidth = 1,linetype = 1 )+
#     geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#     theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
#     NULL)


(nr_t <-  df_res_NR %>%
   slice(range_plot) %>% 
   mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15,
          T_H = H/(gb*0.92*75.38)+Tair,
          d_can_air = Tcan-Tair,
          d_can_leaf = Tcan-tcleaf,
          d_lw_leaf = T_lw-tcleaf,
          d_h_leaf = T_H-tcleaf) %>% 
   pivot_longer(c(d_can_air,d_can_leaf,d_lw_leaf,d_h_leaf)) %>%
   ggplot() +
   geom_abline(slope=0,intercept=0,linetype=2)+
   geom_line(mapping = aes(timestamp,value, color=name),linewidth = 1 )+
   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
   theme_bw()+xlab("Date")+ylab(expression(Delta*"T [ºC]"))+
   scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                 'd_can_leaf' = expression("("*T[can] - T[leaf]*")"),
                                 'd_lw_leaf' = expression("("*T[LW] - T[leaf]*")"),
                                 'd_h_leaf' = expression("("*T[H] - T[leaf]*")")),
                      values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                      name =""
   )+
   theme(legend.text = element_text(size = 16))+
   NULL)

# df_res_NR_filter <- df_res_NR %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res_NR_filter$assim,df_res_NR_filter$GPP)
# df_res2_NR_filter <- df_res2_NR %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res2_NR_filter$assim,df_res2_NR_filter$GPP)

#H
# plot.eval.dens(df_res_NR$Qc,df_res_NR$H)
plot.eval.dens(df_res_NR$tcleaf,(df_res_NR$H/(df_res_NR$gb*0.92*75.38)+df_res_NR$Tair))

#T LW
plot.eval.dens(df_res_NR$tcleaf,(df_res_NR$LW_OUT/(0.96*5.67e-8))^0.25-273.15)

#NETRAD
# plot.eval.dens(df_res_NR$Rnet,df_res_NR$NETRAD)

#T
plot.eval.dens(df_res_NR$tcleaf,df_res_NR$Tcan)

#T
# plot.eval.dens(df_res_NR$Tair,df_res_NR$Tcan)

#USTAR
# plot.eval.dens(df_res_NR$ust,df_res_NR$Ustar)

#LE
# df_res_NR_filter <- df_res_NR %>% filter(SW_IN>25,LAI>2,LE>0)
# plot.eval.dens(df_res_NR_filter$e*18.01528*2230*1e-6,df_res_NR_filter$LE)
# df_res2_NR_filter <- df_res2_NR %>% filter(SW_IN>25,LAI>2,LE>0)
# plot.eval.dens(df_res2_NR_filter$e*18.01528*2230*1e-6,df_res2_NR_filter$LE)
















########################################################################
# 02.Read the data US-Wrc
########################################################################
filename <- filenames.fluxnet[29]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]

data_flx_pre <- fluxnet_data %>% cbind(site = site) %>% mutate(timestamp = lubridate::ymd_hms(DateTime))

fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)

#### read metadata ####
sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

data_flx <- data_flx_pre %>% 
  left_join(sites_metadata) %>% 
  mutate(CO2 = case_when(is.na(CO2)~400,
                         TRUE~CO2)) %>%
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan, Ustar) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP))
summary(data_flx)

res_WR <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                           co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                           ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA,# data_flx$Ustar, 
                           canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                           elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001, 
                           netrad = NA, #data_flx$NETRAD,
                           beta = 146.0, c_cost = 0.41, 
                           do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                           do_acclimation = TRUE, 
                           upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                           acclim_days = 15, weighted_accl = TRUE,
                           energy_params = list(
                             epsleaf = 0.96, #thermal absorptivity of the leaf
                             ste_bolz = 5.67e-8, #W m^-2 K^-4
                             cpm = 75.38, #J mol^-1 ºC-1
                             J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                             frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                             fanir = 0.35 #Fraction of NIR absorbed
                           ))

df_res_WR <- as_tibble(res_WR) %>% cbind(data_flx)

res2_WR <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                            co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                            ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                            canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                            elv = data_flx$elev, z = data_flx$z, leafwidth = 0.001,
                            netrad = NA,#data_flx$NETRAD, 
                            beta = 146.0, c_cost = 0.41,
                            do_leaftemp = FALSE,  gb_method = "Su_2001",
                            do_acclimation = TRUE, 
                            upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                            acclim_days = 15, weighted_accl = TRUE,
                            energy_params = list(
                              epsleaf = 0.98, #thermal absorptivity of the leaf
                              ste_bolz = 5.67e-8, #W m^-2 K^-4
                              cpm = 75.38, #J mol^-1 ºC-1
                              J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                              frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                              fanir = 0.35 #Fraction of NIR absorbed
                            ))

df_res2_WR <- as_tibble(res2_WR) %>% cbind(data_flx)



range_plot <- c(4521:4598)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_WR$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_WR$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

# (wr_no <- ggplot() +
#    geom_line(data = df_res2_WR %>%
#                slice(range_plot),
#              mapping = aes(timestamp,GPP_DT), color = "black")+
#    geom_line(data = df_res2_WR %>%
#                slice(range_plot),
#              mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#    geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#    theme_bw()+xlab("Date")+ylab("GPP")+
#    NULL)
# 
# (wr_yes <- ggplot() +
#     
#     geom_line(data = df_res2_WR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,GPP_DT), color = "black")+
#     geom_line(data = df_res_WR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
#     geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#     theme_bw()+xlab("Date")+ylab("GPP")+
#     NULL)
# 
# (wr_t <- ggplot() +
#     geom_line(data = df_res_WR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
#     geom_line(data = df_res_WR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
#     geom_ribbon(data = df_res_WR %>%
#                   slice(range_plot), 
#                 mapping= aes(x= timestamp , ymin=Tcan-Tcan_sd,ymax=Tcan+Tcan_sd),
#                 color= "grey60",alpha=0.5)+
#     geom_line(data = df_res_WR %>%
#                 slice(range_plot),
#               mapping = aes(timestamp,Tcan), color = "green2",linewidth = 1 )+
#     geom_line(data = df_res_WR %>%
#                 slice(range_plot) %>% 
#                 mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
#               mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
#     geom_line(data = df_res_WR%>%
#                 slice(range_plot) %>% 
#                 mutate(T_H = H/(gb*0.92*75.38)+Tair),
#               mapping = aes(timestamp,T_H), color = "grey40",linewidth = 1,linetype = 1 )+
#     geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#     theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
#     NULL)

(wr_t <-  df_res_WR %>%
   slice(range_plot) %>% 
   mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15,
          T_H = H/(gb*0.92*75.38)+Tair,
          d_can_air = Tcan-Tair,
          d_can_leaf = Tcan-tcleaf,
          d_lw_leaf = T_lw-tcleaf,
          d_h_leaf = T_H-tcleaf) %>% 
   pivot_longer(c(d_can_air,d_can_leaf,d_lw_leaf,d_h_leaf)) %>%
   ggplot() +
   geom_abline(slope=0,intercept=0,linetype=2)+
   geom_line(mapping = aes(timestamp,value, color=name),linewidth = 1 )+
   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
   theme_bw()+xlab("Date")+ylab(expression(Delta*"T [ºC]"))+
   scale_color_manual(labels = c('d_can_air' = expression("("*T[can] - T[air]*")"),
                                   'd_can_leaf' = expression("("*T[can] - T[leaf]*")"),
                                   'd_lw_leaf' = expression("("*T[LW] - T[leaf]*")"),
                                   'd_h_leaf' = expression("("*T[H] - T[leaf]*")")),
                        values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                      name =""
                        )+
   theme(legend.text = element_text(size = 16))+
   NULL)

# df_res_WR_filter <- df_res_WR %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res_WR_filter$assim,df_res_WR_filter$GPP)
# df_res2_WR_filter <- df_res2_WR %>% filter(SW_IN>25,LAI>2)
# plot.eval.dens(df_res2_WR_filter$assim,df_res2_WR_filter$GPP)

#H
# plot.eval.dens(df_res_WR$Qc,df_res_WR$H)
plot.eval.dens(df_res_WR$tcleaf,(df_res_WR$H/(df_res_WR$gb*0.92*75.38)+df_res_WR$Tair))

#T LW
plot.eval.dens(df_res_WR$tcleaf,(df_res_WR$LW_OUT/(0.96*5.67e-8))^0.25-273.15)


#NETRAD
# plot.eval.dens(df_res_WR$Rnet,df_res_WR$NETRAD)

#T
plot.eval.dens(df_res_WR$tcleaf,df_res_WR$Tcan)

#T
# plot.eval.dens(df_res_WR$Tair,df_res_WR$Tcan)

#USTAR
# plot.eval.dens(df_res_WR$ust,df_res_WR$Ustar)

#LE
# df_res_WR_filter <- df_res_WR %>% filter(SW_IN>25)
# plot.eval.dens(df_res_WR_filter$e*18.01528*2230*1e-6,df_res_WR_filter$LE)
# df_res2_WR_filter <- df_res2_WR %>% filter(SW_IN>25)
# plot.eval.dens(df_res2_WR_filter$e*18.01528*2230*1e-6,df_res2_WR_filter$LE)










library(ggpubr)

ggarrange(ha_no,ha_yes,me_no,me_yes,nr_no,nr_yes,wr_no,wr_yes,
          labels = c("US-Ha (without leaf temperature)",
                     "US-Ha (with leaf temperature)   ",
                     "US-Me2 (without leaf temperature)",
                     "US-Me2 (with leaf temperature)   ",
                     "US-NR1 (without leaf temperature)",
                     "US-NR1 (with leaf temperature)   ",
                     "US-Wrc (without leaf temperature)",
                     "US-Wrc (with leaf temperature)   "
          ),
          ncol=1,
          label.x = 0.57,
          align = "v"
)



ggarrange(ha_t,me_t,nr_t,wr_t,
          labels = c("US-Ha",
                     "US-Me2",
                     "US-NR1",
                     "US-Wrc"
          ),
          ncol=1,
          label.x = 0.85,
          align = "v",
          common.legend = TRUE
)






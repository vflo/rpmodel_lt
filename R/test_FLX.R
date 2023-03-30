#### TEST SCRIPT ####

########################################################################
# 01.load the libraries
########################################################################
library(xts)
# library(Evapotranspiration)
# library(readFlux)
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
filenames.fluxnet<- list.files(path="R/data/FLX", "*.csv$", full.names=TRUE,recursive = FALSE)


########################################################################
# 02.Read the data BE-Vie
########################################################################
filename <- filenames.fluxnet[6]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
site
df <- tibble(timestamp=ind) %>% 
  cbind(as_tibble(fluxnet_data)) %>% 
  mutate(date = lubridate::as_date(timestamp))

fapar_df <- read_csv(filenames.fluxnet[1]) %>% 
  dplyr::select(date,spline) %>% 
  dplyr::rename(FAPAR=spline)

lai_df <- read_csv(filenames.fluxnet[11]) %>% 
  mutate(lai = case_when(QA>80~NA_real_,
                         TRUE~LAI)) %>% 
  mutate(time_foo=seq(1,n())) %>%
  mutate(LAI=spline(time_foo,lai ,n=n())$y,
         LAI = case_when(LAI<0~0,
                         TRUE~LAI),
         date = lubridate::as_date(time)) %>%
  dplyr::select(-time_foo,-time)
  
df <- df %>% left_join(fapar_df) %>% left_join(lai_df)


data_flx <- df %>%
  # left_join(sites_metadata) %>%
  mutate(PPFD = PPFD_IN,
         CO2 = CO2_F_MDS,
         Tair = TA_F_MDS,
         VPD = VPD_F_MDS,
         WS = WS_F,
         FAPAR = FAPAR,
         LAI = 2,
         PA = PA_F,
         Tcan = TA_F_MDS,
         elev = 493,
         z = 30,
         Hc = 28) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA)
# summary(data_flx)
data_flx <- data_flx %>% slice_tail(n=10000)

res_BE <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
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
                          epsleaf = 0.96, #thermal absorptivity of the leaf
                          ste_bolz = 5.67e-8, #W m^-2 K^-4
                          cpm = 75.38, #J mol^-1 ºC-1
                          J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                          frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                          fanir = 0.35 #Fraction of NIR absorbed
                        ))

df_res_BE <- as_tibble(res_BE) %>% cbind(data_flx)

res2_BE <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                         co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                         ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                         canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                         elv = data_flx$elev, z = data_flx$z, leafwidth = 0.01,
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

df_res2_BE <- as_tibble(res2_BE) %>% cbind(data_flx)

# res3_BE <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
#                      co2 = data_flx$CO2,fapar = data_flx$FAPAR,
#                      ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
#                      elv = data_flx$elev,do_ftemp_kphio = FALSE
#                     ) %>% as.data.frame()
# 
# df_res3_BE <- as_tibble(res3_BE) %>% cbind(data_flx)


range_plot <- c(4508:4843)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_BE$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_BE$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

be_no <- ggplot() +
  geom_ribbon(data = df_res2_BE %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_BE %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res2_BE %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

be_yes <- ggplot() +
  geom_ribbon(data = df_res2_BE %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_BE %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res_BE %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  # geom_line(data = df_res3_BE %>%
  #             slice(range_plot),
  #           mapping = aes(timestamp,assim), color = "brown")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

be_t <- ggplot() +
  geom_line(data = df_res_BE %>%
              slice(range_plot),
            mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
  geom_line(data = df_res_BE %>%
              slice(range_plot),
            mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
  geom_line(data = df_res_BE %>%
              slice(range_plot) %>% 
              mutate(T_H = H_F_MDS/(gb*0.92*75.38)+Tair),
            mapping = aes(timestamp,T_H), color = "grey40",linewidth = 0.5,linetype = 1)+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
  NULL

plot.eval.dens(df_res_BE$assim,df_res_BE$GPP_DT_CUT_REF)
plot.eval.dens(df_res2_BE$assim,df_res2_BE$GPP_DT_CUT_REF)

#H
plot.eval.dens(df_res_BE$Qc, df_res_BE$H_F_MDS)

#NETRAD
plot.eval.dens(df_res_BE$Rnet,df_res_BE$NETRAD)

#LE
df_res_BE_filter <- df_res_BE %>% filter(SW_IN_F>25)
plot.eval.dens(df_res_BE_filter$e*18.01528*2230*1e-6,df_res_BE_filter$LE_F_MDS)
df_res2_BE_filter <- df_res2_BE %>% filter(SW_IN_F>25)
plot.eval.dens(df_res2_BE_filter$e*18.01528*2230*1e-6,df_res2_BE_filter$LE_F_MDS)






########################################################################
# 02.Read the data CH-Cha
########################################################################
filename <- filenames.fluxnet[7]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
site
df <- tibble(timestamp=ind) %>% 
  cbind(as_tibble(fluxnet_data)) %>% 
  mutate(date = lubridate::as_date(timestamp))

fapar_df <- read_csv(filenames.fluxnet[2]) %>% 
  dplyr::select(date,spline) %>% 
  dplyr::rename(FAPAR=spline)

lai_df <- read_csv(filenames.fluxnet[12]) %>% 
  mutate(lai = case_when(QA>80~NA_real_,
                         TRUE~LAI)) %>% 
  mutate(time_foo=seq(1,n())) %>%
  mutate(LAI=spline(time_foo,lai ,n=n())$y,
         LAI = case_when(LAI<0~0,
                         TRUE~LAI),
         date = lubridate::as_date(time)) %>%
  dplyr::select(-time_foo,-time)

df <- df %>% left_join(fapar_df) %>% left_join(lai_df)


data_flx <- df %>%
  # left_join(sites_metadata) %>%
  mutate(PPFD = PPFD_IN,
         CO2 = CO2_F_MDS,
         Tair = TA_F_MDS,
         VPD = VPD_F_MDS,
         WS = WS_F,
         FAPAR = FAPAR,
         LAI = LAI,
         PA = PA_F,
         Tcan = TA_F_MDS,
         elev = 393,
         z = 2,
         Hc = 1) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA)
# summary(data_flx)
data_flx <- data_flx %>% slice_tail(n=10000)

res_CH <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
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

df_res_CH <- as_tibble(res_CH) %>% cbind(data_flx)

res2_CH <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                            co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                            ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                            canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                            elv = data_flx$elev, z = data_flx$z, leafwidth = 0.01,
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

df_res2_CH <- as_tibble(res2_CH) %>% cbind(data_flx)

# res3_CH <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
#                         co2 = data_flx$CO2,fapar = data_flx$FAPAR,
#                         ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
#                         elv = data_flx$elev,do_ftemp_kphio = FALSE) %>% as.data.frame()
# 
# df_res3_CH <- as_tibble(res3_CH) %>% cbind(data_flx)


range_plot <- c(5014:5350)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_CH$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_CH$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

ch_no <- ggplot() +
  geom_ribbon(data = df_res2_CH %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_CH %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res2_CH %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  # geom_line(data = df_res3_CH %>%
  #             slice(range_plot),
  #           mapping = aes(timestamp,assim), color = "brown")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

ch_yes <- ggplot() +
  geom_ribbon(data = df_res2_CH %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_CH %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res_CH %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

ch_t <- ggplot() +
  geom_line(data = df_res_CH %>%
              slice(range_plot),
            mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
  geom_line(data = df_res_CH %>%
              slice(range_plot),
            mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
  geom_line(data = df_res_CH %>%
              slice(range_plot) %>% 
              mutate(T_H = H_F_MDS/(gb*0.92*75.38)+Tair),
            mapping = aes(timestamp,T_H), color = "grey40",linewidth = 0.5,linetype = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
  NULL

plot.eval.dens(df_res_CH$assim,df_res_CH$GPP_DT_CUT_REF)
plot.eval.dens(df_res2_CH$assim,df_res2_CH$GPP_DT_CUT_REF)

#H
plot.eval.dens(df_res_CH$Qc,df_res_CH$H_F_MDS)

#NETRAD


#LE
df_res_CH_filter <- df_res_CH %>% filter(SW_IN_F>25)
plot.eval.dens(df_res_CH_filter$e*18.01528*2230*1e-6,df_res_CH_filter$LE_F_MDS)
df_res2_CH_filter <- df_res2_CH %>% filter(SW_IN_F>25)
plot.eval.dens(df_res2_CH_filter$e*18.01528*2230*1e-6,df_res2_CH_filter$LE_F_MDS)
















########################################################################
# 02.Read the data FI-Hyy
########################################################################
filename <- filenames.fluxnet[8]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
site
df <- tibble(timestamp=ind) %>% 
  cbind(as_tibble(fluxnet_data)) %>% 
  mutate(date = lubridate::as_date(timestamp))

fapar_df <- read_csv(filenames.fluxnet[3]) %>% 
  dplyr::select(date,spline) %>% 
  dplyr::rename(FAPAR=spline)

lai_df <- read_csv(filenames.fluxnet[13]) %>% 
  mutate(lai = case_when(QA>80~NA_real_,
                         TRUE~LAI)) %>% 
  mutate(time_foo=seq(1,n())) %>%
  mutate(LAI=spline(time_foo,lai ,n=n())$y,
         LAI = case_when(LAI<0~0,
                         TRUE~LAI),
         date = lubridate::as_date(time)) %>%
  dplyr::select(date,LAI)

foo <- fapar_df%>% left_join(lai_df)
df <- df  %>% left_join(foo)


data_flx <- df %>%
  # left_join(sites_metadata) %>%
  mutate(PPFD = PPFD_IN,
         CO2 = CO2_F_MDS,
         Tair = TA_F_MDS,
         VPD = VPD_F,
         WS = WS_F,
         FAPAR = FAPAR,
         LAI = LAI,
         PA = PA_F,
         Tcan = TA_F_MDS,
         elev = 181,
         z = 30,
         Hc = 28) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA)
# summary(data_flx)
data_flx <- data_flx %>% slice_tail(n=10000)

res_FI <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
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

df_res_FI <- as_tibble(res_FI) %>% cbind(data_flx)

res2_FI <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
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

df_res2_FI <- as_tibble(res2_FI) %>% cbind(data_flx)

# res3_FI <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
#                         co2 = data_flx$CO2,fapar = data_flx$FAPAR,
#                         ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
#                         elv = data_flx$elev,do_ftemp_kphio = FALSE) %>% as.data.frame()
# 
# df_res3_FI <- as_tibble(res3_FI) %>% cbind(data_flx)


range_plot <- c(5054:5390)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_FI$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_FI$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

fi_no <- ggplot() +
  geom_ribbon(data = df_res2_FI %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_FI %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res2_FI %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

fi_yes <- ggplot() +
  geom_ribbon(data = df_res2_FI %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_FI %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res_FI %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

fi_t <- ggplot() +
  geom_line(data = df_res_FI %>%
              slice(range_plot),
            mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
  geom_line(data = df_res_FI %>%
              slice(range_plot),
            mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
  geom_line(data = df_res_FI %>%
              slice(range_plot) %>% 
              mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
            mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
  geom_line(data = df_res_FI %>%
              slice(range_plot) %>% 
              mutate(T_H = H_F_MDS/(gb*0.92*75.38)+Tair),
            mapping = aes(timestamp,T_H), color = "grey40",linewidth = 0.5,linetype = 1)+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
  NULL

plot.eval.dens(df_res_FI$assim,df_res_FI$GPP_DT_CUT_REF)
plot.eval.dens(df_res2_FI$assim,df_res2_FI$GPP_DT_CUT_REF)

#H
plot.eval.dens(df_res_FI$Qc,df_res_FI$H_F_MDS)

#NETRAD
plot.eval.dens(df_res_FI$Rnet,df_res_FI$NETRAD)

#LE
df_res_FI_filter <- df_res_FI %>% filter(SW_IN_F>25)
plot.eval.dens(df_res_FI_filter$e*18.01528*2230*1e-6,df_res_FI_filter$LE_F_MDS)
df_res2_FI_filter <- df_res2_FI %>% filter(SW_IN_F>25)
plot.eval.dens(df_res2_FI_filter$e*18.01528*2230*1e-6,df_res2_FI_filter$LE_F_MDS)













########################################################################
# 02.Read the data GF-Guy
########################################################################
filename <- filenames.fluxnet[9]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
site
df <- tibble(timestamp=ind) %>% 
  cbind(as_tibble(fluxnet_data)) %>% 
  mutate(date = lubridate::as_date(timestamp))

fapar_df <- read_csv(filenames.fluxnet[4]) %>% 
  dplyr::select(date,spline) %>% 
  dplyr::rename(FAPAR=spline)

lai_df <- read_csv(filenames.fluxnet[14]) %>% 
  mutate(lai = case_when(QA>80~NA_real_,
                         TRUE~LAI)) %>% 
  mutate(time_foo=seq(1,n())) %>%
  mutate(#LAI=spline(time_foo,lai ,n=n())$y,
         LAI = case_when(LAI<0~0,
                         TRUE~LAI),
         date = lubridate::as_date(time)) %>%
  dplyr::select(-time_foo,-time)

df <- df %>% left_join(fapar_df) %>% left_join(lai_df)


data_flx <- df %>%
  # left_join(sites_metadata) %>%
  mutate(PPFD = PPFD_IN,
         CO2 = CO2_F_MDS,
         Tair = TA_F_MDS,
         VPD = VPD_F_MDS,
         WS = WS_F,
         FAPAR = FAPAR,
         LAI = LAI,
         PA = PA_F,
         Tcan = TA_F_MDS,
         elev = 181,
         z = 30,
         Hc = 28) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA)
# summary(data_flx)
data_flx <- data_flx %>% slice_tail(n=10000)

res_GF <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                           co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                           ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA,# data_flx$Ustar, 
                           canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                           elv = data_flx$elev, z = data_flx$z, leafwidth = 0.3, 
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

df_res_GF <- as_tibble(res_GF) %>% cbind(data_flx)

res2_GF <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
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

df_res2_GF <- as_tibble(res2_GF) %>% cbind(data_flx)
# 
# res3_GF <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
#                         co2 = data_flx$CO2,fapar = data_flx$FAPAR,
#                         ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
#                         elv = data_flx$elev,do_ftemp_kphio = FALSE) %>% as.data.frame()
# 
# df_res3_GF <- as_tibble(res3_GF) %>% cbind(data_flx)


range_plot <- c(5005:5344)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_GF$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_GF$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

gf_no <- ggplot() +
  geom_ribbon(data = df_res2_GF %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_GF %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res2_GF %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

gf_yes <- ggplot() +
  geom_ribbon(data = df_res2_GF %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_GF %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res_GF %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

gf_t <- ggplot() +
  geom_line(data = df_res_GF %>%
              slice(range_plot),
            mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
  geom_line(data = df_res_GF %>%
              slice(range_plot),
            mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
  geom_line(data = df_res_GF %>%
              slice(range_plot) %>% 
              mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
            mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
  geom_line(data = df_res_GF%>%
              slice(range_plot) %>% 
              mutate(T_H = H_F_MDS/(gb*0.92*75.38)+Tair),
            mapping = aes(timestamp,T_H), color = "grey40",linewidth = 0.5,linetype = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
  NULL

plot.eval.dens(df_res_GF$assim,df_res_GF$GPP_DT_CUT_REF)
plot.eval.dens(df_res2_GF$assim,df_res2_GF$GPP_DT_CUT_REF)

#H
plot.eval.dens(df_res_GF$Qc,df_res_GF$H_F_MDS)

#NETRAD
plot.eval.dens(df_res_GF$Rnet,df_res_GF$NETRAD)

#LE
df_res_GF_filter <- df_res_GF %>% filter(SW_IN_F>25)
plot.eval.dens(df_res_GF_filter$e*18.01528*2230*1e-6,df_res_GF_filter$LE_F_MDS)
df_res2_GF_filter <- df_res2_GF %>% filter(SW_IN_F>25)
plot.eval.dens(df_res2_GF_filter$e*18.01528*2230*1e-6,df_res2_GF_filter$LE_F_MDS)










########################################################################
# 02.Read the data US-UMB
########################################################################
filename <- filenames.fluxnet[10]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
site
df <- tibble(timestamp=ind) %>% 
  cbind(as_tibble(fluxnet_data)) %>% 
  mutate(date = lubridate::as_date(timestamp))

fapar_df <- read_csv(filenames.fluxnet[5]) %>% 
  dplyr::select(date,spline) %>% 
  dplyr::rename(FAPAR=spline)

lai_df <- read_csv(filenames.fluxnet[15]) %>% 
  mutate(lai = case_when(QA>80~NA_real_,
                         TRUE~LAI)) %>% 
  mutate(time_foo=seq(1,n())) %>%
  mutate(LAI=spline(time_foo,lai ,n=n(),method = "fmm")$y,
         # LAI=zoo::na.approx(lai,x=time_foo),
         LAI = case_when(LAI<0~0,
                         TRUE~LAI),
         date = lubridate::as_date(time)) %>%
  dplyr::select(-time_foo,-time)

df <- df %>% left_join(fapar_df) %>% left_join(lai_df)

data_flx <- df %>%
  # left_join(sites_metadata) %>%
  mutate(PPFD = PPFD_IN,
         CO2 = CO2_F_MDS,
         Tair = TA_F_MDS,
         VPD = VPD_F_MDS,
         WS = WS_F,
         FAPAR = FAPAR,
         LAI = LAI,
         PA = PA_F,
         Tcan = TA_F_MDS,
         elev = 234,
         z = 46,
         Hc = 20) %>% 
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA) 
# summary(data_flx)
data_flx <- data_flx %>% slice_tail(n=10000)%>% 
  filter(LAI>0)

res_US <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                           co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                           ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA,# data_flx$Ustar, 
                           canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                           elv = data_flx$elev, z = data_flx$z, leafwidth = 0.20, 
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

df_res_US <- as_tibble(res_US) %>% cbind(data_flx)

res2_US <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$Tair, vpd = data_flx$VPD*100,
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

df_res2_US <- as_tibble(res2_US) %>% cbind(data_flx)

# res3_US <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
#                         co2 = data_flx$CO2,fapar = data_flx$FAPAR,
#                         ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
#                         elv = data_flx$elev,do_ftemp_kphio = FALSE) %>% as.data.frame()
# 
# df_res3_US <- as_tibble(res3_US) %>% cbind(data_flx)


range_plot <- c(3372:3514)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res2_US$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res2_US$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

#### GPP plots and analysis ####

us_no <- ggplot() +
  geom_ribbon(data = df_res2_US %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_US %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res2_US %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

us_yes <- ggplot() +
  geom_ribbon(data = df_res2_US %>%
                slice(range_plot), 
              mapping= aes(x= timestamp , ymin=GPP_DT_CUT_05,ymax=GPP_DT_CUT_95),
              color= "grey60",alpha=0.5)+
  geom_line(data = df_res2_US %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP_DT_CUT_REF), color = "black")+
  geom_line(data = df_res_US %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red",linewidth = 1 )+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("GPP")+
  NULL

us_t <- ggplot() +
  geom_line(data = df_res_US %>%
              slice(range_plot),
            mapping = aes(timestamp,Tair), color = "blue",linewidth = 1 )+
  geom_line(data = df_res_US %>%
              slice(range_plot),
            mapping = aes(timestamp,tcleaf), color = "red",linewidth = 1 )+
  geom_line(data = df_res_US %>%
              slice(range_plot) %>% 
              mutate(T_lw = (LW_OUT/(0.96*5.67e-8))^0.25-273.15),
            mapping = aes(timestamp,T_lw), color = "black",linewidth = 1 )+
  geom_line(data = df_res_US%>%
              slice(range_plot) %>% 
              mutate(T_H = H_F_MDS/(gb*0.92*75.38)+Tair),
            mapping = aes(timestamp,T_H), color = "grey40",linewidth = 0.5 ,linetype = 1)+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+xlab("Date")+ylab("Temperature [ºC]")+
  NULL

df_res_US_filter <- df_res_US %>% filter(LAI>0.1)
df_res2_US_filter <- df_res2_US %>% filter(LAI>0.1)
plot.eval.dens(df_res_US_filter$assim,df_res_US_filter$GPP_DT_CUT_REF)
plot.eval.dens(df_res2_US_filter$assim,df_res2_US_filter$GPP_DT_CUT_REF)

#H

plot.eval.dens(df_res_US_filter$Qc,df_res_US_filter$H_F_MDS)

#NETRAD
plot.eval.dens(df_res_US_filter$Rnet,df_res_US_filter$NETRAD)

#LE
df_res_US_filter <- df_res_US %>% filter(SW_IN_F>25,LAI>0.1)
plot.eval.dens(df_res_US_filter$e*18.01528*2230*1e-6,df_res_US_filter$LE_F_MDS)
df_res2_US_filter <- df_res2_US %>% filter(SW_IN_F>25,LAI>0.1)
plot.eval.dens(df_res2_US_filter$e*18.01528*2230*1e-6,df_res2_US_filter$LE_F_MDS)











library(ggpubr)

ggarrange(be_no,be_yes,ch_no,ch_yes,fi_no,fi_yes,gf_no,gf_yes,us_no,us_yes,
          labels = c("BE-Vie (without leaf temperature)",
                     "BE-Vie (with leaf temperature)   ",
                     "CH-cha (without leaf temperature)",
                     "CH-cha (with leaf temperature)   ",
                     "FI-Hyy (without leaf temperature)",
                     "FI-Hyy (with leaf temperature)   ",
                     "GF-Guy (without leaf temperature)",
                     "GF-Guy (with leaf temperature)   ",
                     "US-Umb (without leaf temperature)",
                     "US-Umb (with leaf temperature)   "
                     ),
          ncol=1,
          label.x = 0.57,
          align = "v"
          )



ggarrange(be_t,ch_t,fi_t,gf_t,us_t,
          labels = c("BE-Vie",
                     "CH-cha ",
                     "FI-Hyy",
                     "GF-Guy",
                     "US-Umb"
          ),
          ncol=1,
          label.x = 0.85,
          align = "v"
)












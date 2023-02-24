#### TEST SCRIPT ####

########################################################################
# 01.load the libraries
########################################################################
library(xts)
library(Evapotranspiration)
library(readFlux)
library(ggplot2)
library(zoo)
library(tidyverse)
if(!require(devtools)){install.packages("devtools")}
# devtools::install_github("aukkola/FluxnetLSM", build_vignettes = TRUE)
library("FluxnetLSM")
sapply(list("R/rpmodel_core.R","R/rpmodel.R","R/rpmodel_subdaily.R",
            "R/subroutines.R","R/include_fapar_lai.R","R/include_albedo.R"),source,.GlobalEnv)
library(MASS)
library(ggplot2)
library(viridis)
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
#load the estations md
# stations_md<-read.csv("R/data/fluxnet_list_selection.csv")

########################################################################
# 02.Read the data
########################################################################
# sites<-c('AU-ASM','DE-RuR','ZA-Kru','GH-Ank','US-UMB')

# stations_md<-subset(stations_md,stations_md$SITE_ID %in% sites)
# filename <- filenames.fluxnet[5]
# filename <- filenames.fluxnet[112]
# filename <- filenames.fluxnet[100]
filename <- filenames.fluxnet[31]
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
  na.omit(cols=c("PPFD",'CO2','Tair','FAPAR','VPD','WS','LAI','PA','Tcan')) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(PPFD<0~0,
                         TRUE~GPP))
summary(data_flx)
# data_flx <- data_flx %>% slice(1:4000)

res <- rpmodel_subdaily(TIMESTAMP = data_flx$DateTime, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                        co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                        ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                        elv = data_flx$elev, z = data_flx$z, beta = 146.0, c_cost = 0.41, do_leaftemp = TRUE, do_acclimation = TRUE, 
                        upscaling_method = "max_rad", hour_reference_T = c(10,12,14), acclim_days = 15, weighted_accl = TRUE,
                        energy_params = list(
                          epsleaf = 0.98, #thermal absorptivity of the leaf
                          ste_bolz = 5.67e-8, #W m^-2 K^-4
                          cpm = 75.38, #J mol^-1 ºC-1
                          J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                          frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                          fanir = 0.35 #Fraction of NIR absorbed
                        ))

df_res <- as_tibble(res) %>% cbind(data_flx)

res2 <- rpmodel_subdaily(TIMESTAMP = data_flx$DateTime, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                         co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                        ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                        elv = data_flx$elev, z = data_flx$z, beta = 146.0, c_cost = 0.41, do_leaftemp = FALSE, do_acclimation = TRUE, 
                        upscaling_method = "noon", hour_reference_T = c(10,12,14), acclim_days = 15, weighted_accl = TRUE,
                        energy_params = list(
                          epsleaf = 0.98, #thermal absorptivity of the leaf
                          ste_bolz = 5.67e-8, #W m^-2 K^-4
                          cpm = 75.38, #J mol^-1 ºC-1
                          J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                          frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                          fanir = 0.35 #Fraction of NIR absorbed
                        ))

df_res2 <- as_tibble(res2) %>% cbind(data_flx)

res3 <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
                co2 = data_flx$CO2,fapar = data_flx$FAPAR,
                ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
                elv = data_flx$elev) %>% as.data.frame()

df_res3 <- as_tibble(res3) %>% cbind(data_flx)

res4 <- rpmodel_subdaily(TIMESTAMP = data_flx$DateTime, tc = data_flx$Tcan, vpd = data_flx$VPD_leaf_calc,
                         co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                         ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                         elv = data_flx$elev, z = data_flx$z, beta = 146.0, c_cost = 0.41, do_leaftemp = FALSE, do_acclimation = TRUE, 
                         upscaling_method = "noon", hour_reference_T = c(10,12,14), acclim_days = 15, weighted_accl = TRUE,
                         energy_params = list(
                           epsleaf = 0.98, #thermal absorptivity of the leaf
                           ste_bolz = 5.67e-8, #W m^-2 K^-4
                           cpm = 75.38, #J mol^-1 ºC-1
                           J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                           frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                           fanir = 0.35 #Fraction of NIR absorbed
                         ))

df_res4 <- as_tibble(res4) %>% cbind(data_flx)

range_plot <- c(1090:1200)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res$DateTime) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res$DateTime) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

df_res %>%
  slice(range_plot) %>%
  ggplot() +
  # geom_line(aes(timestamp,T_CANOPY), color = "grey40")+
  # geom_ribbon(aes(timestamp,ymin=T_CANOPY - T_CANOPYsd, ymax=T_CANOPY + T_CANOPYsd),alpha = 0.2)+
  geom_line(aes(DateTime,Tair), color = "red")+
  geom_line(aes(DateTime,tcleaf), color = "blue")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  # geom_col(aes(DateTime,P))+
  # geom_line(aes(DateTime,LEAF_WET/10), color = "purple")+
  ylab("Temperature ºC")+
  theme_bw()+
  NULL
# 
df_res %>%
  # slice(590:700) %>%
  ggplot()+
  geom_point(aes(VPD*100,vpd_leaf))

ggplot() +
  geom_line(data = df_res2 %>%
              slice(range_plot),
            mapping = aes(DateTime,gs), color = "red")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,gs), color = "blue")+
  geom_line(data = df_res3 %>%
              slice(range_plot),
            mapping = aes(DateTime,gs), color = "green")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

#### GPP plots and analysis ####

ggplot() +
  geom_line(data = df_res2 %>%
              slice(range_plot),
            mapping = aes(DateTime,assim), color = "red")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,assim), color = "blue")+
  geom_line(data = df_res3 %>%
              slice(range_plot),
            mapping = aes(DateTime,assim), color = "green")+
  geom_line(data = df_res4 %>%
              slice(range_plot),
            mapping = aes(DateTime,assim), color = "purple")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,GPP), color = "black")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,-GPP_alt), color = "grey40")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

df_res %>% filter(!is.na(GPP)) %>%  mutate(dens = get_density(assim,GPP, n = 100)) %>%
  filter(PPFD>50,P == 0) %>%
  # filter(P == 0) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res2 %>% filter(!is.na(GPP))%>% mutate(dens = get_density(assim,GPP, n = 100)) %>%
  filter(PPFD>50,P == 0) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res4 %>% filter(!is.na(GPP),!is.na(assim))%>% mutate(dens = get_density(assim,GPP, n = 100)) %>%
  filter(PPFD>50,P == 0) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(GPP~assim,
   data = df_res %>% filter(PPFD>50,P == 0)
   ) %>% summary()

lm(GPP~assim,
   data = df_res2 %>% filter(PPFD>50,P == 0)
   ) %>% summary()

lm(GPP~assim,
   data = df_res4 %>% filter(PPFD>50,P == 0)
) %>% summary()

#TLEAF

df_res %>% mutate(
  delta_T = Tair-tcleaf,
  dens = get_density(Tair,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(Tair,tcleaf, color = PPFD))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res %>% mutate(
  delta_T = Tair-tcleaf,
  dens = get_density(Tair,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(Tair,tcleaf, color = e))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res %>% mutate(
  delta_T = Tair-tcleaf,
  dens = get_density(Tair,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(Tair,tcleaf, color = assim))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

# 
# 
# df_res %>% mutate(
#   sigma = 5.67e-8,
#   epsleaf = 0.98,
#   tkleaf = tcleaf + 273.15,
#   tk = Tair + 273.15,
#   Qtirleaf = epsleaf*sigma*tkleaf^4,
#   dens = get_density(Qtirleaf,LW_OUT, n = 100)) %>%
#   filter(PPFD_IN>50,P == 0) %>%
#   ggplot(aes(Qtirleaf,LW_OUT, color =dens))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
#   scale_color_viridis_c(option = "C")
# 
# lm(LW_OUT~Qtirleaf,
#    data = df_res %>% 
#      mutate(sigma = 5.67e-8,
#             epsleaf = 0.98,
#             tkleaf = tcleaf + 273.15,
#             tk = Tair + 273.15,
#             Qtirleaf = epsleaf*sigma*tkleaf^4)%>% filter(PPFD_IN>50,P == 0)
# ) %>% summary()
# 
# #Tair
# df_res %>% 
#   mutate(sigma = 5.67e-8,
#          epsleaf = 0.98,
#          tkleaf = tcleaf + 273.15,
#          tk = Tair + 273.15,
#          Qtirleaf = epsleaf*sigma*tk^4,
#          dens = get_density(Qtirleaf,LW_OUT, n = 100)) %>%
#   filter(PPFD_IN>50,P == 0) %>%
#   ggplot(aes(Qtirleaf,LW_OUT, color =dens))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
#   scale_color_viridis_c(option = "C")
# 
# lm(LW_OUT~Qtirleaf,
#    data = df_res %>% 
#      mutate(sigma = 5.67e-8,
#             epsleaf = 0.98,
#             tkleaf = tcleaf + 273.15,
#             tk = Tair + 273.15,
#             Qtirleaf = epsleaf*sigma*tk^4) %>% filter(PPFD_IN>50,P == 0)
# ) %>% summary()
# 
# 
# ggplot() +
#   geom_line(data = df_res %>%
#               slice(590:800) %>% 
#               mutate(sigma = 5.67e-8,
#                      epsleaf = 0.98,
#                      tkleaf = tcleaf + 273.15,
#                      tk = Tair + 273.15,
#                      Qtirleaf = epsleaf*sigma*tkleaf^4
#                      ),
#             mapping = aes(DateTime,Qtirleaf), color = "blue")+
#   geom_line(data = df_res %>%
#               slice(590:800) %>% 
#               mutate(sigma = 5.67e-8,
#                      epsleaf = 0.98,
#                      tkleaf = tcleaf + 273.15,
#                      tk = Tair + 273.15,
#                      Qtir = epsleaf*sigma*tk^4),
#             mapping = aes(DateTime,Qtir), color = "red")+
#   geom_line(data = df_res %>%
#               slice(590:800),
#             mapping = aes(DateTime,LW_OUT), color = "black")+
#   geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
#   theme_bw()+
#   NULL

# library(psych)
# pairs.panels(data_flx[,-1])

# 
# 
df_res %>% mutate(dens = get_density(Tcan,Tair, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(Tair,Tcan, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

df_res %>% mutate(dens = get_density(Tcan,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(tcleaf,Tcan, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

df_res %>% mutate(dens = get_density(Tair,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(tcleaf,Tair, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

lm(Tcan~tcleaf,data = df_res#%>% filter(PPFD_IN>50,P == 0)
     ) %>% 
  summary()

lm(Tcan~Tair,data = df_res# %>%filter(PPFD_IN>50,P == 0)
   ) %>% summary()


### LE
ggplot() +
  geom_line(data = df_res2 %>%
              slice(range_plot),
            mapping = aes(DateTime,e*18.01528*2230*1e-6), color = "red")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,e*18.01528*2230*1e-6), color = "blue")+
  geom_line(data = df_res3 %>%
              slice(range_plot),
            mapping = aes(DateTime,gs*VPD*100*18.01528*2230*1e-6), color = "green")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,LE), color = "black")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

df_res %>%
  filter(!is.na(e), !is.na(LE)) %>%
  mutate(
  le = e*18.01528*2230*1e-6,
  dens = get_density(LE,le, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(LE,le, color = PPFD))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res2 %>%
  filter(!is.na(e), !is.na(LE)) %>%
  mutate(
    le = e*18.01528*2230*1e-6,
    dens = get_density(LE,le, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(LE,le, color = PPFD))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(LE~le,
   data = df_res %>%
     filter(!is.na(e), !is.na(LE)) %>%
     mutate(le = e*18.01528*2230*1e-6)
) %>% summary()

lm(LE~le,
   data = df_res2 %>%
     filter(!is.na(e), !is.na(LE)) %>%
     mutate(le = e*18.01528*2230*1e-6)
) %>% summary()



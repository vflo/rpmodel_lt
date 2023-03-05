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
# data_flx <- data_flx %>% slice(1000:3000)
#Approximate epsleaf
data_flx %>% 
  filter(SW_IN<25,
         WS<0.01
         ) %>% 
  mutate(epsleaf = LW_OUT/(5.67e-8*(Tcan+273.15)^4)) %>% 
  summarise(epsleaf = mean(epsleaf,na.rm =TRUE))

res <- rpmodel_subdaily(TIMESTAMP = data_flx$DateTime, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                        co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                        ppfd = data_flx$PPFD, u=data_flx$WS, ustar = data_flx$Ustar, 
                        canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                        elv = data_flx$elev, z = data_flx$z, d = 0.001, 
                        netrad = data_flx$NETRAD,
                        beta = 146.0, c_cost = 0.41, 
                        do_leaftemp = TRUE,  gb_method = "Su_2001",
                        do_acclimation = TRUE, 
                        upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                        acclim_days = 15, weighted_accl = TRUE,
                        energy_params = list(
                          epsleaf = 0.9858975, #thermal absorptivity of the leaf
                          ste_bolz = 5.67e-8, #W m^-2 K^-4
                          cpm = 75.38, #J mol^-1 ºC-1
                          J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                          frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                          fanir = 0.35 #Fraction of NIR absorbed
                        ))

df_res <- as_tibble(res) %>% cbind(data_flx) %>% filter(Q_tcleaf == 1)

res2 <- rpmodel_subdaily(TIMESTAMP = data_flx$DateTime, tc = data_flx$Tair, vpd = data_flx$VPD*100,
                         co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                        ppfd = data_flx$PPFD, u=data_flx$WS, ustar = NA, 
                        canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                        elv = data_flx$elev, z = data_flx$z, d = 0.1,
                        netrad = data_flx$NETRAD, beta = 146.0, c_cost = 0.41,
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

df_res2 <- as_tibble(res2) %>% cbind(data_flx) %>% filter(res$Q_tcleaf == 1)

res3 <- rpmodel_core(tc = data_flx$Tair, vpd = data_flx$VPD*100,
                co2 = data_flx$CO2,fapar = data_flx$FAPAR,
                ppfd = data_flx$PPFD, patm =data_flx$PA*1000, 
                elv = data_flx$elev) %>% as.data.frame()

df_res3 <- as_tibble(res3) %>% cbind(data_flx)%>% filter(res$Q_tcleaf == 1)

res4 <- rpmodel_subdaily(TIMESTAMP = data_flx$DateTime, tc = data_flx$Tcan, vpd = data_flx$VPD_leaf_calc,
                         co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                         ppfd = data_flx$PPFD, u=data_flx$WS, ustar = data_flx$Ustar, 
                         canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                         elv = data_flx$elev, z = data_flx$z, d = 0.1, 
                         netrad = data_flx$NETRAD, beta = 146.0, c_cost = 0.41, 
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

df_res4 <- as_tibble(res4) %>% cbind(data_flx)%>% filter(res$Q_tcleaf == 1)

# range_plot <- c(43090:43300)
range_plot <- c(2090:2300)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res$DateTime) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res$DateTime) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

mod <- lm((Tcan-tcleaf)~WS+H+NETRAD+Ustar+LAI+VPD+LE+gb,data = df_res) 
summary(mod)

plot(mod)


df_res %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(LW_OUT,delta_T, n = 100)) %>% 
  ggplot(aes(x=LW_OUT,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(H)) %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(H,delta_T, n = 100)) %>% 
  ggplot(aes(x=H,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(Qc)) %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(Qc,delta_T, n = 100)) %>% 
  ggplot(aes(x=Qc,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(Qc)) %>% 
  mutate(
    delta_T = Tair-tcleaf,
    dens = get_density(Qc,delta_T, n = 100)) %>% 
  ggplot(aes(x=Qc,y=(Tair-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(NETRAD,delta_T, n = 100)) %>% 
  ggplot(aes(x=NETRAD,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")


df_res %>% 
  filter(!is.na(H),!is.na(LE)) %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density((NETRAD - H - LE),delta_T, n = 100)) %>% 
  ggplot(aes(x=(NETRAD - H - LE),y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")


df_res %>% 
  filter(!is.na(H),!is.na(LE)) %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density((NETRAD),Rnet, n = 100)) %>% 
  ggplot(aes(x=(NETRAD),y=(Rnet)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  # filter(!is.na(H)) %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(WS,delta_T, n = 100)) %>% 
  ggplot(aes(x=WS,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(Ustar),SW_IN<25) %>%
  # filter((Tcan-tcleaf)>10) %>% 
  mutate(
    delta_T = Tcan-tcleaf,
    hour = lubridate::hour(timestamp),
    dens = get_density(Ustar,delta_T, n = 100)) %>% 
  ggplot(aes(x=Ustar,y=(Tcan-tcleaf)))+
  geom_point(aes(color=(NETRAD - H - LE)))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(Ustar),SW_IN<25) %>%
  mutate(
    delta_T = Tcan-tcleaf,
    hour = lubridate::hour(timestamp),
    dens = get_density(Ustar,delta_T, n = 100)) %>% 
  ggplot(aes(x=Ustar,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(Ustar),SW_IN>25) %>%
  mutate(
    delta_T = Tair-tcleaf,
    hour = lubridate::hour(timestamp),
    dens = get_density(Ustar,delta_T, n = 100)) %>% 
  ggplot(aes(x=Ustar,y=(Tair-tcleaf)))+
  geom_point(aes(color=NETRAD))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  # filter(!is.na(Ustar)) %>%
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(ust,delta_T, n = 100)) %>% 
  ggplot(aes(x=ust,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  # filter(!is.na(Ustar)) %>%
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(SWC,delta_T, n = 100)) %>% 
  ggplot(aes(x=SWC,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  # filter(!is.na(Ustar)) %>%
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(LAI,delta_T, n = 100)) %>% 
  ggplot(aes(x=LAI,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")

df_res %>% 
  filter(!is.na(LE)) %>%
  mutate(
    delta_T = Tcan-tcleaf,
    dens = get_density(LE,delta_T, n = 100)) %>% 
  ggplot(aes(x=LE,y=(Tcan-tcleaf)))+
  geom_point(aes(color=dens))+
  scale_color_viridis_c(option = "C")






df_res %>%
  slice(range_plot) %>% 
  filter(P==0)%>% 
  ggplot() +
  geom_line(aes(timestamp,Tcan), color = "grey40")+
  # geom_ribbon(aes(timestamp,ymin=T_CANOPY - T_CANOPYsd, ymax=T_CANOPY + T_CANOPYsd),alpha = 0.2)+
  geom_line(aes(DateTime,Tair), color = "red")+
  geom_line(aes(DateTime,tcleaf), color = "blue")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(DateTime,vpd_leaf/100), color = "black")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  # geom_line(aes(DateTime,WS))+
  # geom_line(aes(DateTime,H*0.1))+
  # geom_line(aes(DateTime,LEAF_WET/10), color = "purple")+
  ylab("Temperature ºC")+
  theme_bw()+
  NULL
# 
df_res %>%
  # slice(590:700) %>%
  mutate(
    VPD = VPD*100,
    dens = get_density(VPD,vpd_leaf, n = 100)) %>%
  ggplot(aes(VPD,vpd_leaf, color = dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

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
  # geom_line(data = df_res %>%
  #             slice(range_plot),
  #           mapping = aes(DateTime,-GPP_alt), color = "grey40")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

df_res %>% 
  filter(!is.na(GPP)) %>% 
  filter(PPFD>25) %>%
  mutate(dens = get_density(assim,GPP, n = 100)) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res2 %>% 
  filter(!is.na(GPP))%>%
  filter(PPFD>25) %>%
  mutate(dens = get_density(assim,GPP, n = 100)) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res4 %>% 
  filter(!is.na(GPP),!is.na(assim))%>% 
  filter(PPFD>25) %>%
  mutate(dens = get_density(assim,GPP, n = 100)) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(GPP~assim,
   data = df_res %>% filter(SW_IN>25)
   ) %>% summary()

lm(GPP~assim,
   data = df_res2 %>% filter(SW_IN>25)
   ) %>% summary()

lm(GPP~assim,
   data = df_res4 %>% filter(SW_IN>25)
) %>% summary()

#TLEAF
# 
# df_res %>% mutate(
#   delta_T = Tair-tcleaf,
#   dens = get_density(Tair,tcleaf, n = 100)) %>%
#   # filter(PPFD_IN>50,P == 0) %>%
#   ggplot(aes(Tair,tcleaf, color = PPFD))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
#   scale_color_viridis_c(option = "C")
# 
# df_res %>% mutate(
#   delta_T = Tair-tcleaf,
#   dens = get_density(Tair,tcleaf, n = 100)) %>%
#   # filter(PPFD_IN>50,P == 0) %>%
#   ggplot(aes(Tair,tcleaf, color = e))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
#   scale_color_viridis_c(option = "C")
# 
# df_res %>% mutate(
#   delta_T = Tair-tcleaf,
#   dens = get_density(Tair,tcleaf, n = 100)) %>%
#   # filter(PPFD_IN>50,P == 0) %>%
#   ggplot(aes(Tair,tcleaf, color = assim))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
#   scale_color_viridis_c(option = "C")
# 
# df_res %>% mutate(
#   delta_T = Tair-tcleaf,
#   dens = get_density(Tair,Tcan, n = 100)) %>%
#   filter(PPFD>25) %>%
#   ggplot(aes(Tair,Tcan, color = assim))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)+
#   geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
#   scale_color_viridis_c(option = "C")


df_res %>%  
  # filter(PPFD>25) %>%
  mutate(dens = get_density(Tcan,Tair, n = 100)) %>%
  ggplot(aes(Tair,Tcan, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

df_res %>%
  # filter(PPFD>25) %>% 
  mutate(dens = get_density(Tcan,tcleaf, n = 100)) %>%
  ggplot(aes(tcleaf,Tcan, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

df_res %>%
  # filter(PPFD>25) %>% 
  mutate(dens = get_density(Tair,tcleaf, n = 100)) %>%
  ggplot(aes(tcleaf,Tair, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

lm(Tcan~tcleaf,data = df_res# %>% filter(SW_IN>25)
     ) %>% 
  summary()

lm(Tcan~Tair,data = df_res #%>%filter(SW_IN>25)
   ) %>% summary()


### H
df_res %>%  
  filter(PPFD>25) %>%
  filter(!is.na(H)) %>%
  mutate(dens = get_density(Qc, H, n = 100)) %>%
  ggplot(aes(Qc, H, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

### QTIR
df_res %>%  
  # filter(PPFD>25, P == 0) %>% 
  filter(!is.na(Qtirleaf)) %>%
  mutate(dens = get_density(Qtirleaf, LW_OUT, n = 100)) %>%
  ggplot(aes(Qtirleaf, LW_OUT, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

df_res %>%  
  # filter(PPFD>25) %>%
  filter(!is.na(Qtirleaf)) %>%
  mutate(
    Qtircan = 0.98*5.67e-8*(Tcan+273.15)^4,
    dens = get_density(Qtirleaf, Qtircan, n = 100)) %>%
  ggplot(aes(Qtirleaf, Qtircan, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")


df_res2 %>%  
  # filter(PPFD>25, P == 0) %>% 
  # filter(!is.na(Qtirleaf)) %>%
  mutate(
    Qtirleaf = 0.98*5.67e-8*(Tair+273.15)^4,
    dens = get_density(Qtirleaf, LW_OUT, n = 100)) %>%
  ggplot(aes(Qtirleaf, LW_OUT, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  scale_color_viridis_c(option = "C")

lm(LW_OUT~Qtirleaf,
   data = df_res %>%
     filter(!is.na(e), !is.na(LE)) %>%
     # filter(PPFD>25,P == 0) %>%
     mutate(le = e*18.01528*2230*1e-6)
) %>% summary()

lm(LW_OUT~Qtirleaf,
   data = df_res2 %>%
     # filter(!is.na(e), !is.na(LE)) %>%
     # filter(PPFD>25,P == 0) %>%
     mutate(Qtirleaf = 0.98*5.67e-8*(Tair+273.15)^4,)
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
  filter(SW_IN>25) %>%
  mutate(
  le = e*18.01528*2230*1e-6,
  dens = get_density(LE,le, n = 100)) %>%
  ggplot(aes(le,LE, color = dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res2 %>%
  filter(!is.na(e), !is.na(LE)) %>%
  filter(SW_IN>25) %>%
  mutate(
    le = e*18.01528*2230*1e-6,
    dens = get_density(LE,le, n = 100)) %>%
  ggplot(aes(le,LE, color = dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res4 %>%
  filter(!is.na(e), !is.na(LE)) %>%
  filter(SW_IN>25) %>%
  mutate(
    le = e*18.01528*2230*1e-6,
    dens = get_density(LE,le, n = 100)) %>%
  ggplot(aes(le,LE, color = dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(LE~lE,
   data = df_res %>%
     filter(!is.na(e), !is.na(LE)) %>%
     filter(SW_IN>25) %>%
     mutate(le = e*18.01528*2230*1e-6)
) %>% summary()

lm(LE~le,
   data = df_res2 %>%
     filter(SW_IN>25) %>%
     filter(PPFD>25,P == 0) %>%
     mutate(le = e*18.01528*2230*1e-6)
) %>% summary()

lm(LE~le,
   data = df_res4 %>%
     filter(!is.na(e), !is.na(LE)) %>%
     filter(SW_IN>25) %>%
     mutate(le = e*18.01528*2230*1e-6)
) %>% summary()

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
sapply(list("R/rpmodel_core.R","R/rpmodel.R","R/rpmodel_subdaily.R","R/subroutines.R"),source,.GlobalEnv)

########################################################################
# # Reading Eddy Covariance data using the readFlux Package
########################################################################
########################################################################
# 01.find the filename
########################################################################
#load files available
#### read the files' paths
filenames.fluxnet<- list.files(path="R/data", "*.csv$", full.names=TRUE,recursive = TRUE)
#load the estations md
# stations_md<-read.csv("R/data/fluxnet_list_selection.csv")

########################################################################
# 02.Read the data
########################################################################
# sites<-c('AU-ASM','DE-RuR','ZA-Kru','GH-Ank','US-UMB')

# stations_md<-subset(stations_md,stations_md$SITE_ID %in% sites)
filename <- filenames.fluxnet[2]
fluxnet_data <- data.table::fread(filename, header = T, 
                              quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                              colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M", 
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]

data_flx <- tibble(timestamp=ind) %>% cbind(as_tibble(fluxnet_data))

# data_flx_GH_Ank<-readFluxdata(filename=filenames.fluxnet[4],elev = stations_md$elv[4],0.1)

# CO2_maunaLoa<-read.table('https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt')
# CO2_maunaLoa[CO2_maunaLoa==-9.99]<-NA
# names(CO2_maunaLoa)<-c('year','month','decimal_date','average_CO2', 'interpolated','trend','#days')
# CO2_maunaLoa <- CO2_maunaLoa %>% dplyr::select(year, month, average_CO2)
# soil_site <- getSoilSite(stations_md$LAT[3],stations_md$LON[3])

data_flx <- data_flx %>% 
  mutate(T_CANOPY = rowMeans(select(data_flx,starts_with("T_CANOPY")),na.rm=TRUE),
         TA = rowMeans(select(data_flx,starts_with("TA_")),na.rm=TRUE),
         VPD = rowMeans(select(data_flx,starts_with("VPD")),na.rm=TRUE),
         PPFD_IN = rowMeans(select(data_flx,starts_with("PPFD_IN")),na.rm=TRUE),
         LEAF_WET = rowMeans(select(data_flx,starts_with("LEAF_WET")),na.rm=TRUE),
         WS = rowMeans(select(data_flx,starts_with("WS")),na.rm=TRUE),
         SW_IN = rowMeans(select(data_flx,starts_with("SW_IN")),na.rm=TRUE),
         CO2 = rowMeans(select(data_flx,starts_with("CO2")),na.rm=TRUE),
         PA = rowMeans(select(data_flx,starts_with("PA_")),na.rm=TRUE),
         SWC = SWC_2_3_1,
         LE = LE_2_1_1) %>%
  dplyr::select(timestamp,T_CANOPY,TA,VPD,PPFD_IN,LEAF_WET,WS,SW_IN,SWC,LE,CO2,PA) %>% 
  na.omit()
  # filter(!is.na(T_CANOPY),!is.na(TA),!is.na(VPD),!is.na(PPFD_IN),!is.na)


# library(psych)
# pairs.panels(data_flx[,-1])

# data_flx %>% 
#   ggplot(aes(T_CANOPY,TA))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)
# 
# data_flx %>% 
#   ggplot(aes(LE,TA-T_CANOPY))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)
# 
# # mod <- lm(TA-T_CANOPY~VPD+LE+PPFD_IN+LEAF_WET+WS+SW_IN,data=data_flx)
# # summary(mod)
# # step(mod) %>% summary()
# # relaimpo::calc.relimp(mod)
# # car::vif(mod)
# # plot(mod)
# 
# 
# 
# 
# data_flx %>% filter(SWC>=0) %>% 
#   ggplot()+
#   geom_line(aes(timestamp,SWC))
# 
# data_flx %>% 
#   slice(c(2000:2500)) %>% 
#   ggplot()+
#   geom_line(aes(timestamp,T_CANOPY),color="green")+
#   # geom_line(aes(timestamp,T_CANOPY_2_1_2),color="red")+
#   geom_line(aes(timestamp,TA),color="blue")+
#   geom_line(aes(timestamp,LEAF_WET/10),color="purple")+
#   NULL
# 
# data_flx %>% filter(T_CANOPY>-100 | T_CANOPY>-100,
#                     TA_1_1_1>-100#,SWC_1_1_1>=0
# ) %>% 
#   ggplot(aes((T_CANOPY_2_1_1+T_CANOPY_2_1_2)/2,TA_1_1_1))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)
# 
# data_flx %>% filter(T_CANOPY_1_1_1>-100,TA_1_1_1>-100#,SWC_1_1_1>=0
#                     ) %>% 
#   ggplot()+
#   geom_line(aes(timestamp,LE_1_1_1),color="red")
#########################################################################
# 03.Run pmodel_lt
#########################################################################
# test_fn <- function(x){
#   rpmodel(tc = x$tc, vpd = x$VPD, co2 = x$CO2, fapar = 0.7, 
#           ppfd = x$PPFD/(60*60*12)*1e6, 
#           u= 2,canopy_height = 10, elv = 0, do_leaftemp = TRUE)
# }
# 
# data_flx_GH_Ank %>% 
#   fortify.zoo %>%
#   as_tibble %>% 
#   rename(timestamp = "Index") %>% 
#   filter(!is.na(VPD),!is.na(tc),!is.na(PPFD)) %>% 
#   mutate(year = lubridate::year(timestamp),
#          month = lubridate::month(timestamp)) %>% 
#   left_join(CO2_maunaLoa, by = c("year", "month")) %>% 
#   mutate(CO2 = dplyr::coalesce(CO2, average_CO2)) %>% 
#   slice(1:300) %>%
#   # split(seq(nrow(.))) %>%
#   # rowwise() %>% 
#   mutate(p = purrr::pmap(list(tc, VPD, CO2, PPFD), 
#                               ~rpmodel(
#                                 tc = ..1, vpd = ..2, co2 = ..3, fapar = 0.7, 
#                                 ppfd = ..4/(60*60*24)*1e6, 
#                                 u= NA,canopy_height = 10, elv = 0, do_leaftemp = FALSE
#                               )),
#          lt = purrr::pmap(list(tc, VPD, CO2, PPFD), 
#                          ~rpmodel(
#                            tc = ..1, vpd = ..2, co2 = ..3, fapar = 0.7, 
#                            ppfd = ..4/(60*60*24)*1e6, 
#                            u= NA,canopy_height = 10, elv = 0, do_leaftemp = TRUE
#                          )))%>% 
#   unnest_wider(p)%>% 
#   unnest_wider(lt,names_sep="_")%>%
#   ggplot()+
#   geom_line(aes(timestamp,tc),color="red")+
#   geom_line(aes(timestamp,lt_tcleaf),color="black")+
#   # geom_line(aes(timestamp,gpp),color="red")+
#   # geom_line(aes(timestamp,lt_gpp),color="black")+
#   # geom_abline(intercept=0,slope=1)
#   NULL




res <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$TA, vpd = data_flx$VPD*1000,co2 = data_flx$CO2,fapar = 0.9,
                 ppfd = data_flx$PPFD_IN,u=data_flx$WS,canopy_height = 20, patm =data_flx$PA*1000, elv = 0)





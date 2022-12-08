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
sapply(list("R/rpmodel_core.R","R/rpmodel.R","R/subroutines.R"),source,.GlobalEnv)

  ######################################################################
# # Reading Eddy Covariance data using the readFlux Package
########################################################################
########################################################################
# 01.find the filename
########################################################################
#load files available
#### read the files' paths
filenames.fluxnet<- list.files(path="R/data/", pattern = ".*.FULLSET_H.*.csv$", full.names=TRUE)
#load the estations md
stations_md<-read.csv("R/data/fluxnet_list_selection.csv")

########################################################################
# 02.Read the data
########################################################################
sites<-c('AU-ASM','DE-RuR','ZA-Kru','GH-Ank')

stations_md<-subset(stations_md,stations_md$SITE_ID %in% sites)
data_flx_GH_Ank<-readFluxdata(filename=filenames.fluxnet[4],elev = stations_md$elv[4],0.1)

CO2_maunaLoa<-read.table('https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt')
CO2_maunaLoa[CO2_maunaLoa==-9.99]<-NA
names(CO2_maunaLoa)<-c('year','month','decimal_date','average_CO2', 'interpolated','trend','#days')
CO2_maunaLoa <- CO2_maunaLoa %>% dplyr::select(year, month, average_CO2)
# soil_site <- getSoilSite(stations_md$LAT[3],stations_md$LON[3])


#########################################################################
# 03.Run pmodel_lt
#########################################################################
test_fn <- function(x){
  rpmodel(tc = x$tc, vpd = x$VPD, co2 = x$CO2, fapar = 0.7, 
          ppfd = x$PPFD/(60*60*12)*1e6, 
          u= 2,canopy_height = 10, elv = 0, do_leaftemp = TRUE)
}

data_flx_GH_Ank %>% 
  fortify.zoo %>%
  as_tibble %>% 
  rename(timestamp = "Index") %>% 
  filter(!is.na(VPD),!is.na(tc),!is.na(PPFD)) %>% 
  mutate(year = lubridate::year(timestamp),
         month = lubridate::month(timestamp)) %>% 
  left_join(CO2_maunaLoa, by = c("year", "month")) %>% 
  mutate(CO2 = dplyr::coalesce(CO2, average_CO2)) %>% 
  slice(1:300) %>%
  # split(seq(nrow(.))) %>%
  # rowwise() %>% 
  mutate(p = purrr::pmap(list(tc, VPD, CO2, PPFD), 
                              ~rpmodel(
                                tc = ..1, vpd = ..2, co2 = ..3, fapar = 0.7, 
                                ppfd = ..4/(60*60*24)*1e6, 
                                u= NA,canopy_height = 10, elv = 0, do_leaftemp = FALSE
                              )),
         lt = purrr::pmap(list(tc, VPD, CO2, PPFD), 
                         ~rpmodel(
                           tc = ..1, vpd = ..2, co2 = ..3, fapar = 0.7, 
                           ppfd = ..4/(60*60*24)*1e6, 
                           u= NA,canopy_height = 10, elv = 0, do_leaftemp = TRUE
                         )))%>% 
  unnest_wider(p)%>% 
  unnest_wider(lt,names_sep="_")%>%
  ggplot()+
  geom_line(aes(timestamp,tc),color="red")+
  geom_line(aes(timestamp,lt_tcleaf),color="black")+
  # geom_line(aes(timestamp,gpp),color="red")+
  # geom_line(aes(timestamp,lt_gpp),color="black")+
  # geom_abline(intercept=0,slope=1)
  NULL

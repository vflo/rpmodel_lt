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
filenames.fluxnet<- list.files(path="R/data", "*.csv$", full.names=TRUE,recursive = TRUE)
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
filename <- filenames.fluxnet[124]
fluxnet_data <- data.table::fread(filename, header = T, 
                              quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                              colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M", 
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]

data_flx_pre <- tibble(timestamp=ind) %>% cbind(as_tibble(fluxnet_data)) %>% cbind(si_code = site)

fapar <- read.csv(file="R/data/FAPAR_sites.csv")
fapar_noaa <- read.csv(file="R/data/FAPAR_sites_noaa.csv")
albedo <- read.csv(file="R/data/albedo_sites.csv")

data_flx_pre <- include_fapar_lai(data_flx_pre,fapar_noaa,fapar)
# data_flx_pre <- include_albedo(data_flx_pre,albedo)
# data_flx_GH_Ank<-readFluxdata(filename=filenames.fluxnet[4],elev = stations_md$elv[4],0.1)

# CO2_maunaLoa<-read.table('https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt')
# CO2_maunaLoa[CO2_maunaLoa==-9.99]<-NA
# names(CO2_maunaLoa)<-c('year','month','decimal_date','average_CO2', 'interpolated','trend','#days')
# CO2_maunaLoa <- CO2_maunaLoa %>% dplyr::select(year, month, average_CO2)
# soil_site <- getSoilSite(stations_md$LAT[3],stations_md$LON[3])

# data_flx <- data_flx %>% 
#   mutate(T_CANOPYsd = matrixStats::rowSds(as.matrix(select(data_flx,starts_with("T_CANOPY"))),na.rm=TRUE),
#          T_CANOPY = rowMeans(select(data_flx,starts_with("T_CANOPY")),na.rm=TRUE),
#          TA = rowMeans(select(data_flx,starts_with("TA_")),na.rm=TRUE),
#          VPD = rowMeans(select(data_flx,starts_with("VPD")),na.rm=TRUE),
#          PPFD_IN = rowMeans(select(data_flx,starts_with("PPFD_IN")),na.rm=TRUE),
#          LEAF_WET = rowMeans(select(data_flx,starts_with("LEAF_WET")),na.rm=TRUE),
#          WS = rowMeans(select(data_flx,starts_with("WS")),na.rm=TRUE),
#          SW_IN = rowMeans(select(data_flx,starts_with("SW_IN")),na.rm=TRUE),
#          CO2 = rowMeans(select(data_flx,starts_with("CO2")),na.rm=TRUE),
#          PA = rowMeans(select(data_flx,starts_with("PA_")),na.rm=TRUE),
#          LE = rowMeans(select(data_flx,starts_with("LE_")),na.rm=TRUE),
#          P = P_1_1_1,
#          SWC = SWC_2_3_1) %>%
#   dplyr::select(timestamp,T_CANOPY,T_CANOPYsd,TA,VPD,PPFD_IN,LEAF_WET,WS,SW_IN,SWC,LE,CO2,PA,P) %>% 
#   na.omit(cols=c(1:13)) #%>% slice(595)
  # filter(!is.na(T_CANOPY),!is.na(TA),!is.na(VPD),!is.na(PPFD_IN),!is.na)

#GF-Guy
data_flx <- data_flx_pre %>%
  dplyr::select(-c(paste0("TS_F_MDS_",c(1:2),"_QC"), paste0("SWC_F_MDS_",c(1:4),"_QC")))
#DE-Tha
# data_flx <- data_flx_pre %>%
#   dplyr::select(-c(paste0("TS_F_MDS_",c(1:6),"_QC"), paste0("SWC_F_MDS_",c(1:2),"_QC")))
#CH-Cha
# data_flx <- data_flx_pre %>%
#   dplyr::select(-c(paste0("TS_F_MDS_",c(1:9),"_QC"), paste0("SWC_F_MDS_",c(1:3,5),"_QC")))
data_flx <- data_flx %>% 
  mutate(TA = TA_F_MDS,
         TS_F_MDS = rowMeans(dplyr::select(data_flx,starts_with("TS_")),na.rm=TRUE),
         VPD = VPD_F_MDS,
         PPFD_IN = PPFD_IN,
         LW_IN_F_MDS = LW_IN_F,
         # LW_OUT = LW_OUT,
         USTAR = USTAR,
         H_F_MDS = H_F_MDS,
         WS = WS_F,
         SW_IN = SW_IN_F_MDS,
         CO2 = CO2_F_MDS,
         PA = PA,
         LE = LE_F_MDS,
         GPP = GPP_DT_VUT_MEAN,
         P = P_F,
         SWC = rowMeans(dplyr::select(data_flx,starts_with("SWC_F_")),na.rm=TRUE)) %>%
  dplyr::select(timestamp,TA, TS_F_MDS, VPD, PPFD_IN, LW_IN_F_MDS,# LW_OUT, 
                H_F_MDS, WS, USTAR,
                SW_IN, SWC, LE, CO2, PA, P, GPP, LAI, FAPAR, si_code) %>% 
  na.omit(cols=c(1:19))

sites_metadata <- read.csv(file="R/data/sites_metadata.csv")

data_flx <- data_flx %>% left_join(sites_metadata)

data_flx <- data_flx %>% slice(1:4000)

res <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$TA, vpd = data_flx$VPD*100,
                        co2 = data_flx$CO2, fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                        ppfd = data_flx$PPFD_IN,u=data_flx$WS, ustar = NA, canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                        elv = data_flx$elev, z = data_flx$z, do_leaftemp = TRUE, do_acclimation = TRUE, 
                        upscaling_method = "noon", acclim_days = 15, weighted_accl = TRUE,
                        energy_params = list(
                          epsleaf = 0.98, #thermal absorptivity of the leaf
                          ste_bolz = 5.67e-8, #W m^-2 K^-4
                          cpm = 75.38, #J mol^-1 ºC-1
                          kfFEC = 2.0, #Photon flux to energy μmol J-1 (Meek et al., 1984)
                          fanir = 0.35 #Fraction of NIR absorbed
                        ))

df_res <- as_tibble(res) %>% cbind(data_flx)

res2 <- rpmodel_subdaily(TIMESTAMP = data_flx$timestamp, tc = data_flx$TA, vpd = data_flx$VPD*100,
                         co2 = data_flx$CO2,fapar = data_flx$FAPAR, LAI = data_flx$LAI,
                        ppfd = data_flx$PPFD_IN,u=data_flx$WS, ustar = NA, canopy_height = data_flx$Hc, patm =data_flx$PA*1000, 
                        elv = data_flx$elev, z = data_flx$z, do_leaftemp = FALSE, do_acclimation = TRUE, 
                        upscaling_method = "noon", acclim_days = 15, weighted_accl = TRUE,
                        energy_params = list(
                          epsleaf = 0.98, #thermal absorptivity of the leaf
                          ste_bolz = 5.67e-8, #W m^-2 K^-4
                          cpm = 75.38, #J mol^-1 ºC-1
                          kfFEC = 2.0, #Photon flux to energy μmol J-1 (Meek et al., 1984)
                          fanir = 0.35 #Fraction of NIR absorbed
                        ))

df_res2 <- as_tibble(res2) %>% cbind(data_flx)

res3 <- rpmodel_core(tc = data_flx$TA, vpd = data_flx$VPD*100,
                co2 = data_flx$CO2,fapar = data_flx$FAPAR,
                ppfd = data_flx$PPFD_IN, patm =data_flx$PA*1000, 
                elv = data_flx$elev) %>% as.data.frame()

df_res3 <- as_tibble(res3) %>% cbind(data_flx)

range_plot <- c(1090:1200)
dateRanges <- data.frame(
  start = lubridate::as_date(df_res$timestamp) %>% lubridate::as_datetime() + 18*60*60,
  end   = lubridate::as_date(df_res$timestamp) %>% lubridate::as_datetime() + 30*60*60 )%>%
  slice(range_plot)%>% 
  dplyr::distinct()

df_res %>%
  slice(range_plot) %>%
  ggplot() +
  # geom_line(aes(timestamp,T_CANOPY), color = "grey40")+
  # geom_ribbon(aes(timestamp,ymin=T_CANOPY - T_CANOPYsd, ymax=T_CANOPY + T_CANOPYsd),alpha = 0.2)+
  geom_line(aes(timestamp,TA), color = "red")+
  geom_line(aes(timestamp,tcleaf), color = "blue")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  # geom_col(aes(timestamp,P))+
  # geom_line(aes(timestamp,LEAF_WET/10), color = "purple")+
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
            mapping = aes(timestamp,gs), color = "red")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(timestamp,gs), color = "blue")+
  geom_line(data = df_res3 %>%
              slice(range_plot),
            mapping = aes(timestamp,gs), color = "green")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

#### GPP plots and analysis ####

ggplot() +
  geom_line(data = df_res2 %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "red")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "blue")+
  geom_line(data = df_res3 %>%
              slice(range_plot),
            mapping = aes(timestamp,assim), color = "green")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(timestamp,GPP), color = "black")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

df_res %>% mutate(dens = get_density(assim,GPP, n = 100)) %>%
  # filter(PPFD_IN>10,P == 0) %>%
  # filter(P == 0) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res2 %>% mutate(dens = get_density(assim,GPP, n = 100)) %>%
  # filter(PPFD_IN>10,P == 0) %>%
  ggplot(aes(assim,GPP, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(GPP~assim,
   data = df_res #%>% filter(PPFD_IN>10,P == 0)
   ) %>% summary()

lm(GPP~assim,
   data = df_res2 #%>% filter(PPFD_IN>10,P == 0)
   ) %>% summary()


#TLEAF

df_res %>% mutate(
  delta_T = TA-tcleaf,
  dens = get_density(TA,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(TA,tcleaf, color = PPFD_IN))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res %>% mutate(
  delta_T = TA-tcleaf,
  dens = get_density(TA,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(TA,tcleaf, color = e))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

df_res %>% mutate(
  delta_T = TA-tcleaf,
  dens = get_density(TA,tcleaf, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(TA,tcleaf, color = assim))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")



df_res %>% mutate(
  sigma = 5.67e-8,
  epsleaf = 0.98,
  tkleaf = tcleaf + 273.15,
  tk = TA + 273.15,
  Qtirleaf = epsleaf*sigma*tkleaf^4,
  dens = get_density(Qtirleaf,LW_OUT, n = 100)) %>%
  filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(Qtirleaf,LW_OUT, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(LW_OUT~Qtirleaf,
   data = df_res %>% 
     mutate(sigma = 5.67e-8,
            epsleaf = 0.98,
            tkleaf = tcleaf + 273.15,
            tk = TA + 273.15,
            Qtirleaf = epsleaf*sigma*tkleaf^4)%>% filter(PPFD_IN>50,P == 0)
) %>% summary()

#TA
df_res %>% 
  mutate(sigma = 5.67e-8,
         epsleaf = 0.98,
         tkleaf = tcleaf + 273.15,
         tk = TA + 273.15,
         Qtirleaf = epsleaf*sigma*tk^4,
         dens = get_density(Qtirleaf,LW_OUT, n = 100)) %>%
  filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(Qtirleaf,LW_OUT, color =dens))+
  geom_point()+
  geom_abline(intercept=0,slope=1)+
  geom_smooth(method = "lm", color = "grey20", linetype = 2,alpha=0.1)+
  scale_color_viridis_c(option = "C")

lm(LW_OUT~Qtirleaf,
   data = df_res %>% 
     mutate(sigma = 5.67e-8,
            epsleaf = 0.98,
            tkleaf = tcleaf + 273.15,
            tk = TA + 273.15,
            Qtirleaf = epsleaf*sigma*tk^4) %>% filter(PPFD_IN>50,P == 0)
) %>% summary()


ggplot() +
  geom_line(data = df_res %>%
              slice(590:800) %>% 
              mutate(sigma = 5.67e-8,
                     epsleaf = 0.98,
                     tkleaf = tcleaf + 273.15,
                     tk = TA + 273.15,
                     Qtirleaf = epsleaf*sigma*tkleaf^4
                     ),
            mapping = aes(timestamp,Qtirleaf), color = "blue")+
  geom_line(data = df_res %>%
              slice(590:800) %>% 
              mutate(sigma = 5.67e-8,
                     epsleaf = 0.98,
                     tkleaf = tcleaf + 273.15,
                     tk = TA + 273.15,
                     Qtir = epsleaf*sigma*tk^4),
            mapping = aes(timestamp,Qtir), color = "red")+
  geom_line(data = df_res %>%
              slice(590:800),
            mapping = aes(timestamp,LW_OUT), color = "black")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

# library(psych)
# pairs.panels(data_flx[,-1])

# 
# 
# df_res %>% mutate(dens = get_density(T_CANOPY,TA, n = 100)) %>%
#   filter(PPFD_IN>50,P == 0) %>%  
#   ggplot(aes(TA,T_CANOPY, color =dens))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)
# 
# df_res %>% mutate(dens = get_density(T_CANOPY,tcleaf, n = 100)) %>%
#   filter(PPFD_IN>50,P == 0) %>% 
#   ggplot(aes(tcleaf,T_CANOPY, color =dens))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)
# 
# df_res %>% mutate(dens = get_density(TA,tcleaf, n = 100)) %>%
#   filter(PPFD_IN>50,P == 0) %>% 
#   ggplot(aes(tcleaf,TA, color =dens))+
#   geom_point()+
#   geom_abline(intercept=0,slope=1)
# 
# lm(T_CANOPY~tcleaf,data = df_res%>% 
#      filter(PPFD_IN>50,P == 0)) %>% summary()
# 
# lm(T_CANOPY~TA,data = df_res %>% 
#      filter(PPFD_IN>50,P == 0)) %>% summary()


### LE
ggplot() +
  geom_line(data = df_res2 %>%
              slice(range_plot),
            mapping = aes(timestamp,e*18.01528*2230*1e-6), color = "red")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(timestamp,e*18.01528*2230*1e-6), color = "blue")+
  geom_line(data = df_res3 %>%
              slice(range_plot),
            mapping = aes(timestamp,gs*VPD*100*18.01528*2230*1e-6), color = "green")+
  geom_line(data = df_res %>%
              slice(range_plot),
            mapping = aes(timestamp,LE), color = "black")+
  geom_rect(data = dateRanges, mapping = aes(xmin= start , xmax=end,ymin=-Inf,ymax=Inf), alpha=0.1, fill="black")+
  theme_bw()+
  NULL

df_res %>%
  filter(!is.na(e), !is.na(LE)) %>%
  mutate(
  le = e*18.01528*2230*1e-6,
  dens = get_density(LE,le, n = 100)) %>%
  # filter(PPFD_IN>50,P == 0) %>%
  ggplot(aes(LE,le, color = PPFD_IN))+
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
  ggplot(aes(LE,le, color = PPFD_IN))+
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


##### FUNCTIONS FAPAR LAI #######
library(readr)
library(sf)
library(sp)
library(rgdal)
library(tidyverse)
library(data.table)
library(R.utils)
sites_coord <- read_csv("R/data/sites_coord.csv")
sites_coord <- sites_coord %>% 
  rename(lon = si_long,
         lat = si_lat) %>% 
  dplyr::select(-si_code)

# Define coordinate reference system
prj4string <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
my.projection <- st_crs(prj4string)

# Create sf object
lat_long_sf <- st_as_sf(sites_coord, coords = c("lon", "lat"), crs = my.projection)
st_crs(lat_long_sf)

plot(lat_long_sf)

# Export shapefile
st_write(lat_long_sf, "R/data/sites_coord/sites_coord.shp", driver="ESRI Shapefile")


# modification of dist_merge to retain coordinates of y
dist <- function(x1, y1, x2, y2) {
  ((x1-x2)^2 + (y1-y2)^2)^0.5
}

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


coord <-read_csv("R/data/sites_coord.csv")

#### MODIS ####
fapar <-read_csv("R/data/FAPAR.csv")

fapar %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo

FAPAR_reverse<-dist_merge_mod(coord, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

FAPAR <- foo %>% 
  left_join(FAPAR_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))

write_csv(FAPAR, file="R/data/FAPAR_sites.csv")


#### NOAA ####
fapar <-read_csv("R/data/FAPAR_noaa.csv")

fapar %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat),
         QC = intToBin(FparLai_QC),
         FparLai_QC = str_sub(QC,-2,-1),
         FparLai_QC = strtoi(FparLai_QC, base = 2))->foo

FAPAR_reverse<-dist_merge_mod(coord, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

FAPAR <- foo %>% 
  dplyr::select(-QC) %>% 
  left_join(FAPAR_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long))

write_csv(FAPAR, file="R/data/FAPAR_sites_noaa.csv")

#### ALBEDO ####
albedo <-read_csv("R/data/albedo.csv")

albedo %>% 
  mutate(coord = sub("{\"type\":\"Point\",\"coordinates\":[","",`.geo`,fixed = TRUE),
         coord = sub("]}","",coord)) %>% 
  separate(coord,c("si_long",'si_lat'),",") %>% 
  dplyr::select(-`.geo`) %>% 
  distinct() %>% 
  mutate(si_long = as.numeric(si_long),
         si_lat = as.numeric(si_lat))->foo

albedo_reverse<-dist_merge_mod(coord, foo, 'si_long', 'si_lat', 'si_long', 'si_lat')

albedo <- foo %>% 
  left_join(albedo_reverse%>% 
              dplyr::select(si_code, y_lat, y_long) %>% 
              rename(si_lat = y_lat,
                     si_long = y_long)) %>% 
  filter(BRDF_Albedo_Band_Mandatory_Quality_vis == 0)

write_csv(albedo, file="R/data/albedo_sites.csv")



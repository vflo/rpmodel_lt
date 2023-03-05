size_max_temp_max <- function(d,df,max_temp){
  print(paste("Leaf width:",d))
  
  df <- df %>%
    mutate(month = lubridate::month(DateTime)) %>% 
    group_by(month) %>% 
    mutate(T_month = mean(Tair,na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(T_month_max = max(T_month,na.rm = TRUE)) %>% 
    filter(T_month == T_month_max)
  # q90 <- quantile(df$Tair, probs = 0.95)
  # foo <- subset(df, Tair >= q90, select=timestamp)
  # dates <- foo %>% mutate(dates = lubridate::as_date(timestamp))
  # dates <- dates %>% mutate(dates_less =dates-lubridate::days(15))
  # index <- dates %>% 
  #   split(seq(nrow(.))) %>%
  #   purrr::map(function(x){
  #     tibble(seq_time = seq(from=lubridate::as_datetime(x$dates_less), by=30*60, to=lubridate::as_datetime(x$dates)))
  #     
  #   }) %>% bind_rows() %>% unique()
  # 
  # # df <- df %>% mutate(date = lubridate::as_date(DateTime))
  # df <- subset(df, Tair > q90)

  res <- rpmodel_subdaily(TIMESTAMP = df$DateTime, tc = df$Tair, vpd = df$VPD*100,
                          co2 = df$CO2, fapar = df$FAPAR, LAI = df$LAI,
                          ppfd = df$PPFD, u=df$WS, ustar = NA, 
                          canopy_height = df$Hc, patm =df$PA*1000, 
                          elv = df$elev, z = df$z, d = d, 
                          netrad = df$NETRAD, beta = 146.0, c_cost = 0.41, 
                          do_leaftemp = TRUE, simple_gb = TRUE, do_acclimation = TRUE, 
                          upscaling_method = "max_rad", hour_reference_T = c(10,12,14), 
                          acclim_days = 15, weighted_accl = TRUE,
                          energy_params = list(
                            epsleaf = 0.98, #thermal absorptivity of the leaf
                            ste_bolz = 5.67e-8, #W m^-2 K^-4
                            cpm = 75.38, #J mol^-1 ºC-1
                            J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                            frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                            fanir = 0.35 #Fraction of NIR absorbed
                          ))
  max_x <- max(res$tcleaf,na.rm=TRUE)
  y <- (max_x - max_temp)^2
  print(paste("T_max:",max_x," / RMSE:", y))
  y
  
}



size_max_temp_min <- function(d,df,min_temp){
  print(paste("Leaf width:",d))
  df <- df %>%
    mutate(month = lubridate::month(DateTime)) %>% 
    group_by(month) %>% 
    mutate(T_month = mean(Tair,na.rm = TRUE)) %>%
    filter(T_month>=5)
  
  res <- rpmodel_subdaily(TIMESTAMP = df$DateTime, tc = df$Tair, vpd = df$VPD*100,
                          co2 = df$CO2, fapar = df$FAPAR, LAI = df$LAI,
                          ppfd = df$PPFD, u=df$WS, ustar = NA, 
                          canopy_height = df$Hc, patm =df$PA*1000, 
                          elv = df$elev, z = df$z, d = d, 
                          netrad = df$NETRAD, beta = 146.0, c_cost = 0.41, 
                          do_leaftemp = TRUE, simple_gb = TRUE, do_acclimation = TRUE, 
                          upscaling_method = "max_rad", hour_reference_T = c(10,12,14), 
                          acclim_days = 15, weighted_accl = TRUE,
                          energy_params = list(
                            epsleaf = 0.98, #thermal absorptivity of the leaf
                            ste_bolz = 5.67e-8, #W m^-2 K^-4
                            cpm = 75.38, #J mol^-1 ºC-1
                            J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                            frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                            fanir = 0.35 #Fraction of NIR absorbed
                          ))
  min_x <- min(res$tcleaf,na.rm = TRUE)
  y <- (min_x - min_temp)^2
  print(paste("T_min:",min_x," / RMSE:", y))
  y
  
}











########################################################################
# 01.load the libraries
########################################################################
library(xts)
# library(Evapotranspiration)
# library(readFlux)
library(ggplot2)
library(zoo)
library(tidyverse)
# if(!require(devtools)){install.packages("devtools")}
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
filename <- filenames.fluxnet[1]
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
  drop_na(PPFD,CO2,Tair,FAPAR,VPD,WS,LAI,PA,Tcan) %>% 
  mutate(VPD_leaf_calc = calc_new_vpd(Tcan,Tair,VPD*100),
         GPP = case_when(GPP<0~0,
                         TRUE~GPP))
summary(data_flx)
data_flx <- data_flx %>% slice(1:20000)







opt_max <- optimise(f=size_max_temp_max,
         interval=c(0.0005,3),
         df=data_flx,
         max_temp=50)


opt_min <- optimise(f=size_max_temp_min,
                    interval=c(0.0005,0.05),
                    df=data_flx,
                    min_temp=-5)

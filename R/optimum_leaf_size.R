#### LEAF SIZE OPTIMIZATION ####

########################################################################
# 01.load the libraries
########################################################################
# setwd('/rds/general/user/vflosier/home/opt_leaf')
# getwd()
setwd("/data/r_users/vflo/rpmodel_lt")
library(bigleaf)
library(zoo)
library(tidyverse)
sapply(list("R/rpmodel_core.R","R/rpmodel.R","R/rpmodel_subdaily.R",
            "R/subroutines.R","R/include_fapar_lai.R","R/include_albedo.R"),source,.GlobalEnv)



########################################################################
# 02. load data
########################################################################


#### read the files' paths
filenames <- list.files(path="R/data/df_opt", "*.csv$", full.names=TRUE,recursive = TRUE)

opt_fn <- function(par, df){
  leaf_width <- par
  res <- rpmodel_subdaily(TIMESTAMP = df$timestamp, tc = df$Tair, vpd = df$vpd,
                          co2 = df$co2, fapar = df$FAPAR, LAI = df$LAI,
                          ppfd = df$ppfd, u=df$ws, ustar = NA,# df$Ustar, 
                          canopy_height = df$veg_height, patm =df$pa, 
                          elv = df$elev, z = df$veg_height+2, leafwidth = leaf_width, 
                          netrad = NA, #df$NETRAD,
                          beta = 146.0, c_cost = 0.41, 
                          do_leaftemp = TRUE,  gb_method = "simple",#Choudhury_1988",
                          do_acclimation = TRUE, 
                          upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                          acclim_days = 15, weighted_accl = TRUE, epsleaf = df$Emis_31, #thermal absorptivity of the leaf
                          energy_params = list(
                            ste_bolz = 5.67e-8, #W m^-2 K^-4
                            cpm = 75.38, #J mol^-1 ºC-1
                            J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                            frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                            fanir = 0.35 #Fraction of NIR absorbed
                          ))
  anet <- res$assim - res$rd
  rd <- res$rd
  a <- res$assim
  print(paste("Site:",unique(df$Site),"/ leaf width:",
              leaf_width,"/ Net assimilation:", mean(anet,na.rm=TRUE),
              "/ Assimilation:", mean(a,na.rm=TRUE),
              "/ Rd:", mean(rd,na.rm=TRUE)))
  mean(anet,na.rm=TRUE)
  
}


calc_opt_leaf_size <- function(file_site){
  
  df_site <- read_csv(file=file_site)
  df_site <- df_site %>% 
    dplyr::rowwise() %>% 
    mutate(timestamp=lubridate::with_tz(timestamp,time_zone),
           timestamp = lubridate::force_tz(timestamp,"GMT")) %>% 
    ungroup()
  df_site <- df_site %>% 
    mutate(year=lubridate::year(Date)) %>%
    filter(!is.na(FAPAR),!is.na(LAI), as.character(year) == "2010")
  opt_2010 <- optim(
    par = 0.05,
    fn = opt_fn,
    lower = 0.0001,
    upper = 2,
    method = "L-BFGS-B",
    control = list(fnscale= -1),
    df = df_site
  )
  
  df_site <- read_csv(file=file_site)
  df_site <- df_site %>% 
    dplyr::rowwise() %>% 
    mutate(timestamp=lubridate::with_tz(timestamp,time_zone),
           timestamp = lubridate::force_tz(timestamp,"GMT")) %>% 
    ungroup()
  df_site <- df_site %>% 
    mutate(year=lubridate::year(Date)) %>%
    filter(!is.na(FAPAR),!is.na(LAI), as.character(year) == "2011")
  opt_2011 <- optim(
    par = 0.05,
    fn = opt_fn,
    lower = 0.0001,
    upper = 2,
    method = "L-BFGS-B",
    control = list(fnscale= -1),
    df = df_site
  )
  
  df_site <- read_csv(file=file_site)
  df_site <- df_site %>% 
    dplyr::rowwise() %>% 
    mutate(timestamp=lubridate::with_tz(timestamp,time_zone),
           timestamp = lubridate::force_tz(timestamp,"GMT")) %>% 
    ungroup()
  df_site <- df_site %>% 
    mutate(year=lubridate::year(Date)) %>%
    filter(!is.na(FAPAR),!is.na(LAI), as.character(year) == "2012")
  opt_2012 <- optim(
    par = 0.05,
    fn = opt_fn,
    lower = 0.0001,
    upper = 2,
    method = "L-BFGS-B",
    control = list(fnscale= -1),
    df = df_site
  )
  
  df_site <- read_csv(file=file_site)
  df_site <- df_site %>% 
    dplyr::rowwise() %>% 
    mutate(timestamp=lubridate::with_tz(timestamp,time_zone),
           timestamp = lubridate::force_tz(timestamp,"GMT")) %>% 
    ungroup()
  df_site <- df_site %>% 
    mutate(year=lubridate::year(Date)) %>%
    filter(!is.na(FAPAR),!is.na(LAI), as.character(year) == "2013")
  opt_2013 <- optim(
    par = 0.05,
    fn = opt_fn,
    lower = 0.0001,
    upper = 2,
    method = "L-BFGS-B",
    control = list(fnscale= -1),
    df = df_site
  )
  
  df_site <- read_csv(file=file_site)
  df_site <- df_site %>% 
    dplyr::rowwise() %>% 
    mutate(timestamp=lubridate::with_tz(timestamp,time_zone),
           timestamp = lubridate::force_tz(timestamp,"GMT")) %>% 
    ungroup()
  df_site <- df_site %>% 
    mutate(year=lubridate::year(Date)) %>%
    filter(!is.na(FAPAR),!is.na(LAI), as.character(year) == "2014")
  opt_2014 <- optim(
    par = 0.05,
    fn = opt_fn,
    lower = 0.0001,
    upper = 2,
    method = "L-BFGS-B",
    control = list(fnscale= -1),
    df = df_site
  )
  
  df_site <- read_csv(file=file_site)
  df_site <- df_site %>% 
    dplyr::rowwise() %>% 
    mutate(timestamp=lubridate::with_tz(timestamp,time_zone),
           timestamp = lubridate::force_tz(timestamp,"GMT")) %>% 
    ungroup()
  df_site <- df_site %>% 
    mutate(year=lubridate::year(Date)) %>%
    filter(!is.na(FAPAR),!is.na(LAI), as.character(year) == "2015")
  opt_2015 <- optim(
    par = 0.05,
    fn = opt_fn,
    lower = 0.0001,
    upper = 2,
    method = "L-BFGS-B",
    control = list(fnscale= -1),
    df = df_site
  )
  
  res_df_2010 <- tibble(Site = gsub("/","_",unique(df_site$Site)),
                        leaf_width = opt_2010$par,
                        anet = opt_2010$value,
                        year = 2010)
  
  res_df_2011 <- tibble(Site = gsub("/","_",unique(df_site$Site)),
                        leaf_width = opt_2011$par,
                        anet = opt_2011$value,
                        year = 2011)
  
  res_df_2012 <- tibble(Site = gsub("/","_",unique(df_site$Site)),
                        leaf_width = opt_2012$par,
                        anet = opt_2012$value,
                        year = 2012)
  
  res_df_2013 <- tibble(Site = gsub("/","_",unique(df_site$Site)),
                        leaf_width = opt_2013$par,
                        anet = opt_2013$value,
                        year = 2013)
  
  res_df_2014 <- tibble(Site = gsub("/","_",unique(df_site$Site)),
                        leaf_width = opt_2014$par,
                        anet = opt_2014$value,
                        year = 2014)
  
  res_df_2015 <- tibble(Site = gsub("/","_",unique(df_site$Site)),
                        leaf_width = opt_2015$par,
                        anet = opt_2015$value,
                        year = 2015)
  
  res_df <- bind_rows(res_df_2010,res_df_2011,
                      res_df_2012,res_df_2013,
                      res_df_2014,res_df_2015)

  write.csv(res_df,file=paste0('R/data/results_01/',unique(res_df$Site),".csv"))
  return(res_df)
}


calc_opt_leaf_size(filenames[547])
# envir_vars<-as.character(ls())


# library(doSNOW)
# cl <- makeCluster(2,'SOCK')
# registerDoSNOW(cl) 
# snow::clusterEvalQ(cl,lapply(c('bigleaf','zoo','tidyverse'), library, character.only = TRUE))
# snow::clusterExport(cl, list=envir_vars,envir=environment())
# obs_tsfs<-snow::clusterMap(cl = cl, fun=calc_opt_leaf_size,file_site=filenames[c(73,)])	
# save(obs_tsfs,file='R/data/results_opt_leaf_01.RData')
# stopCluster(cl)




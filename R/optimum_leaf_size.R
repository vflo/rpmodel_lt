#### LEAF SIZE OPTIMIZATION ####

########################################################################
# 01.load the libraries
########################################################################
library(bigleaf)
library(zoo)
library(tidyverse)
sapply(list("R/rpmodel_core.R","R/rpmodel.R","R/rpmodel_subdaily.R",
            "R/subroutines.R","R/include_fapar_lai.R","R/include_albedo.R"),source,.GlobalEnv)



########################################################################
# 02. load data
########################################################################

#### read the files' paths
filenames <- list.files(path="R/data/final_sites", "*.csv$", full.names=TRUE,recursive = TRUE)

opt_fn <- function(leaf_width, df){
  
  res <- rpmodel_subdaily(TIMESTAMP = df$timestamp, tc = df$Tair, vpd = df$VPD*100,
                          co2 = df$CO2, fapar = df$FAPAR, LAI = df$LAI,
                          ppfd = df$PPFD, u=df$WS, ustar = NA,# df$Ustar, 
                          canopy_height = df$Hc, patm =df$PA*1000, 
                          elv = df$elev, z = df$z, leafwidth = leaf_width, 
                          netrad = NA, #df$NETRAD,
                          beta = 146.0, c_cost = 0.41, 
                          do_leaftemp = TRUE,  gb_method = "Choudhury_1988",
                          do_acclimation = TRUE, 
                          upscaling_method = "noon", hour_reference_T = c(10,12,14), 
                          acclim_days = 15, weighted_accl = TRUE,
                          energy_params = list(
                            epsleaf = df$epsleaf, #thermal absorptivity of the leaf
                            ste_bolz = 5.67e-8, #W m^-2 K^-4
                            cpm = 75.38, #J mol^-1 ºC-1
                            J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
                            frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
                            fanir = 0.35 #Fraction of NIR absorbed
                          ))
  
  mean(res$assim,na.rm=TRUE)
  
}


calc_opt_leaf_size <- function(file_site){
  
  df_site <- read.csv(file=file_site)
  opt <- optimr::optimr(
    par = 0.05,
    fn = opt_fn,
    df = df_site,
    lower = 0.0001,
    upper = 2
  )
  
  res_df <- tibble(leaf_width = opt$par,)
  write.csv()
  
}



foreach::foreach(file_site = iter(filenames, by='row')) %dopar% 
  calc_opt_leaf_size(file_site)




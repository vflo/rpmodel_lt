# This code is used to calculate leaf area index based on environmental factors
# Author: Shengchao Qiao, Tsinghua University
# Created date: 24th August 2020
# last updated date: 24th August 2020
##----------------------------------------------
source('PC_model.R') # load the functions used to calculate GPP
library('lamW')


##########-------------------------------------------------------------------------
# function 1: calculate LAI_max based on environmental factors
##########-------------------------------------------------------------------------
cal_LAI_c <- function(Ta,PPFD_ac,VPD,CO2,elv=NA,patm=NA,alpha, index=NA){
  #-----------------------------------------------------------------------
  ## Input:
  # Ta: mean daily air temperature over growing season, degree C
  # PPFD_ac: total PPFD over growing season, mol photon/m2
  # VPD: mean daily vapour pressure deficit, kPa
  # CO2: annual CO2 concentration, ppm
  # elv: elevation, mm
  # patm: mean daily air pressure, kPa
  # alpha: moisture index, the ratio of AET to EET, dimensionless
  # index: the preset maximum value of LAI, 
  
  ## Output:
  # LAI_c: maximum LAI during growing season controled by carbon balance
  
  ##Features: calculate the maximum value of LAI during growing season based on carbon mass-balance
  #-----------------------------------------------------------------------
  
  ## define constant
  LMA <- 35.7 # leaf mass area, g dry mass/m2
  k <- 0.5 # canopy light extinction coefficient
  f_tb <- 0.5 # the ratio of total biomass to GPP
  f_leaf <-0.5 # the ratio of leaf biomass to aboveground biomass before leaf senescent
  f_alpha_g <- -0.23 # the sensitivity coefficient of root ratio (root biomass/ total biomass) to alpha
  
  # the mean value of root ratio to total biomass before lead senescent
  f_root_g <- f_alpha_g*alpha+0.4
  # the fraction of leaf to accumulated GPP 
  eta <- f_tb*(1-f_root_g)*f_leaf
  # calculate LAI
  LUE <- cal_lue(Ta=Ta,VPD = VPD,CO2 = CO2,elv = elv,patm = patm)
  mu <- 2.5*0.5*LUE*PPFD_ac*eta/LMA
  LAI_ma <- mu+2*lambertW0(-k*mu*exp(-k*mu))
  if (is.na(index)){
    LAI_c <- LAI_ma
  } else {
    LAI_c <- min(LAI_ma,runif(1,(index-0.1),(index+0.1)))
  }
  
  
  return(LAI_max)
}

##########-------------------------------------------------------------------------
cal_LAI_w <- function(Ta,PPFD_ac,VPD,CO2,elv=NA,patm=NA,pre, index=0.5){
  #-----------------------------------------------------------------------
  ## Input:
  # Ta: mean daily air temperature over growing season, degree C
  # PPFD_ac: total PPFD over growing season, mol photon/m2
  # VPD: mean daily vapour pressure deficit, kPa
  # CO2: annual CO2 concentration, ppm
  # elv: elevation, mm
  # patm: mean daily air pressure, kPa
  # pre: total precipitation over growing season, mm
  # index: the fraction of average of fAPAR over growing season to maximum fAPAR，default value is 0.5 
  
  ## Output:
  # LAI_w: maximum LAI during growing season controlled by water balance
  
  ##Features: calculate the maximum value of LAI during growing season based on water mass-balance
  #-----------------------------------------------------------------------
  
  ## define constant
  Rue <- 0.9 # the fraction of evapotranspiration (ET) to precipitation
  f_T_to_ET <- 0.7 # the fraction of transpiration (T) to evapotranspiration (ET)
  
  # calculate ET based on precipitation
  ET <- Rue*pre
  T <- ET*f_T_to_ET
  
  # calculate assimilation based on Fick's law
  Ca <- cal_co2_to_ca(co2 = CO2,patm = patm*1000)
  gstar_Ca <- cal_gstar_gepisat(temp = Ta)/Ca
  f1 <- exp(-0.0227*(Ta-25))
  xi <-sqrt(146*(cal_k(temp = Ta,patm = patm*1000)+cal_gstar_gepisat(temp = Ta))/(1.6*f1))
  X <- gstar_Ca+(1-gstar_Ca)*xi/(xi+sqrt(VPD*1000))
  A <- 12*T*1000*Ca*(1-X)/(1.6*VPD*1000*18)
  
  # calculate LAI based on Beer's law
  LUE <- cal_lue(Ta=Ta,VPD = VPD,CO2=CO2,patm = patm)
  fAPAR <- A/(PPFD_ac*LUE)/index
  fAPAR <- ifelse(fAPAR>0.99,0.99,fAPAR)
  LAI_w <- log(1-fAPAR)/(-0.5)
  return(LAI_w)
}


##########-------------------------------------------------------------------------
# function 2: calculate moving average in certain step size
##########-------------------------------------------------------------------------
cal_moving <- function(data,mo=NA){
  #-----------------------------------------------------------------------
  ## Input:
  # data: original data need to calculate moving average
  # mo: the step size of moving average
  
  ## Output:
  # data_mo: the new data by moving average
  
  ##Features: the functiton is used to calculate the moving average in certain step size
  #-----------------------------------------------------------------------
  
    n_data <- length(data) # the number of dataset
    data_mo <- matrix(data = NA,nrow=n_data,ncol = 1)
    for (n in 1:n_data) {
      if(n<mo){
        data_mo[n,1]<-mean(data[1:n])
      } else {
        data_mo[n,1]<-mean(data[(n-mo+1):n])
      }
    }
    
  return(data_mo)
}

##########-------------------------------------------------------------------------
# function 3: calculate the phenological scalar (fPHU) based on air temperature(fPHU)
##########-------------------------------------------------------------------------
cal_fPHU <- function(Ta,Tmin=NA,Tmax=NA,ty=NA){
  #-----------------------------------------------------------------------
  ## Input:
  # Ta: air temperature, degree C
  # Tmin: the base tempetature above which wheat starts grow, degree C
  # Tmax: the upper limit of temperature accumulation in one day, degree C
  # ty: wheat type, spring wheat (sw), winter wheat (ww)
  
  ## Output:
  # PHU: the total heat units over growing season (from planting to harvest), degree C day
  # fPHU: the fraction of accumulated heat units to the total heat units, dimensionless
  
  ## Features: this function is used to calculate the phenological scalar based on air temperature
  #  the fPHU represents the phenological development process of wheat
  #-----------------------------------------------------------------------
  # define constant
  fLAI_1 <- 0.05 # the fraction of LAImax corresponding to 1st point on the optimal LAI development curve 
  fLAI_2 <- 0.95 # the fraction of LAImax corresponding to 2nd point on the optimal LAI development curve 
  
  n <- length(Ta)
  data <- data.frame(Ta_or = Ta)
  
  for (i in 1:n) {
    if (Ta[i]<=Tmin){
      TT <- 0
    } else if (Ta[i]<=26){
      TT <- Ta[i]
    } else if (Ta[i]<=Tmax){
      TT <- (26/8)*(Tmax-Ta[i])
    } else {
      TT <- 0
    }
    
    if (i==1){
      data$phu[i] <- TT
    } else {
      data$phu[i] <- data$phu[i-1]+TT
    }
    
  }
  data$fphu <- data$phu/data$phu[n]

  # different leaf area development parameters for different wheat type.
  # These parameters are from Table A-4 in User's manual of SWAT version2000
  if (ty=='ww'){
    fphu_1 <- 0.05
    fphu_2 <- 0.45
    fphu_sen <- 0.5
  } else {
    fphu_1 <- 0.15
    fphu_2 <- 0.5
    fphu_sen <- 0.6
  }
  # calculate the two inflection point of LAI growing curve during the growth phase
  l2 <- (log(fphu_1/fLAI_1 - fphu_1)-log(fphu_2/fLAI_2 - fphu_2))/(fphu_2-fphu_1)
  l1 <- log(fphu_1/fLAI_1 - fphu_1)+l2*fphu_1
  
  d_sen <- min(which(data$fphu>=fphu_sen))-1
  
  for (i in 1:n) {
    
    if (data$fphu[i]<=fphu_sen ){
      data$fLAI[i]<-data$fphu[i]/(data$fphu[i]+exp(l1-l2*data$fphu[i]))
    } else {
      data$fLAI[i]<-(1-data$fphu[i])^2/(1-fphu_sen)^2*(1-0)
    }
    
  }
  result <- list(phenology=data,type=ty,period=c(ver_end,d_sen))
  return(result)
}

























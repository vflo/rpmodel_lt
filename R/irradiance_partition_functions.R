irradiance_partition <- function(sun_shade = c("sun","shade"), TIMESTAMP, PPFD, LAI, PA, 
                                 PA0 = 101325, fa = 0.426, canopy_reflexion = 0.15, 
                                 rho_cd = 0.036, kd_prime = 0.719, lat, long){
  
  list.of.packages <- c("lubridate", "solrad")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  stopifnot("TIMESTAMP must be a POSIXct, POSIXtl or a Date object" = all(lubridate::is.instant(TIMESTAMP)))
  if(any(PA <10000)) warning("Check if PA is in Pascal")
  if(any(PPFD > 3000)) warning("Check if PPFD is in micromol m-2 s-1")
  
  #Solar elevation angle
  DOY <- lubridate::yday(TIMESTAMP)
  hour <- lubridate::hour(TIMESTAMP)
  DOY_dec <- solrad::DayOfYear(TIMESTAMP)
  beta_angle <- solrad::Altitude(DOY = DOY_dec, 
                                 Lat = lat,
                                 Lon = long,
                                 SLon = long,
                                 DS = 0) *pi/180
  # beam extintion coeficient
  kb = 0.5/sin(beta_angle)
  kb = ifelse(beta_angle<=0, 1e-10,kb )
  # beam and scattered beam extintion coeficient
  kb_prime = 0.46/sin(beta_angle)
  kb_prime = ifelse(beta_angle<=0, 1e-10,kb_prime )
  #fraction of diffuse radiation
  m = (PA/PA0)/sin(beta_angle)
  fd = (1-0.72^m)/(1+0.72^m*(1/fa - 1))
  #beam irradiance horizontal leaves
  rho_h = (1-(1-canopy_reflexion)^0.5)/(1+(1-canopy_reflexion)^0.5)
  #beam irradiance, uniform leaf-angle distribution
  rho_cb = (1-exp(-2*rho_h*kb/(1+kb)))
  #diffuse irradiance
  I_d = PPFD*fd
  I_d = ifelse(I_d<0,0,I_d)
  #beam irradiance
  I_b = PPFD*(1-fd)
  #scattered beam irradiance
  I_bs = I_b*((1-rho_cb)*kb_prime*exp(-kb_prime*LAI)-(1-canopy_reflexion)*kb*exp(-kb*LAI))
  #Irradiance sun exposed
  I_c = (1-rho_cb)*I_b*(1-exp(-kb_prime*LAI))+(1-rho_cd)*I_d*(1-exp(-kd_prime*LAI))
  
  a = I_b*(1-canopy_reflexion)*(1-exp(-kb*LAI))
  b = I_d*(1-rho_cd)*(1-exp(-(kd_prime+kb)*LAI))*kd_prime/(kd_prime+kb)
  c = I_b*(((1-rho_cb)*(1-exp(-(kb_prime+kb)*LAI))*(kb_prime/(kb_prime+kb)))-
            (1-canopy_reflexion)*((1-exp(-2*kb*LAI)))/2)
  
  Isun = a+b+c
  
  Ishade = I_c-Isun
  
  if(sun_shade == "sun"){return(Isun)}else
  if(sun_shade == "shade"){return(Ishade)}else{"sun_shade variable must be either sun or shade"}
  
}

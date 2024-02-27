# ## this document includes the fundamental functions of PCmodel to calculate light use efficiency(LUE),
# #  gross primary productivity (GPP), aboveground biomass (AB) and grain yield (yield).
# ##
# ## updated date: 28th August 2020
# ## by Shengchao Qiao (qsc17@tsinghua.edu.cn) and Han WANG
# 
# ##===================================
# ## define functions
# ##===================================
# 
# # $1. calculate air pressure in Pa
# cal_patm <- function( elv ){
#   #-----------------------------------------------------------------------
#   # Input:    - elevation, m (elv)
#   # Output:   - float, atmospheric pressure at elevation 'elv', Pa (patm)
#   # Features: Returns the atmospheric pressure as a function of elevation
#   #           and standard atmosphere (1013.25 hPa)
#   # Depends:  - connect_sql
#   #           - flux_to_grid
#   #           - get_data_point
#   #           - get_msvidx
#   # Ref:      Allen et al. (1998)
#   #-----------------------------------------------------------------------
#   
#   # Define constants:
#   kPo <- 101325   # standard atmosphere, Pa (Allen, 1973)
#   kTo <- 298.15   # base temperature, K (Prentice, unpublished)
#   kL <- 0.0065    # temperature lapse rate, K/m (Allen, 1973)
#   kG <- 9.80665   # gravitational acceleration, m/s^2 (Allen, 1973)
#   kR <- 8.3143    # universal gas constant, J/mol/K (Allen, 1973)
#   kMa <- 0.028963 # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
#   
#   # Convert elevation to pressure, Pa:
#   patm <- kPo*(1.0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
#   
#   return (patm)
# }
# 
# # $2. calculate K (MM coefficient of Rubisco) in Pa
# cal_k <- function(temp, patm) {
#   #-----------------------------------------------------------------------
#   # Input:    - float, air temperature, deg C (temp)
#   #           - float, atmospheric pressure, Pa (patm)
#   # Output:   float, Pa (mmk)
#   # Features: Returns the temperature & pressure dependent Michaelis-Menten
#   #           coefficient, K (Pa).
#   # Ref:      Bernacchi et al. (2001), Improved temperature response 
#   #           functions for models of Rubisco-limited photosynthesis, 
#   #           Plant, Cell and Environment, 24, 253--259.
#   #-----------------------------------------------------------------------
#   
#   # Define constants
#   kc25 <- 39.97      # Pa, assuming 25 deg C & 98.716 kPa
#   ko25 <- 2.748e4    # Pa, assuming 25 deg C & 98.716 kPa
#   dhac <- 79430      # J/mol
#   dhao <- 36380      # J/mol
#   kR   <- 8.3145     # J/mol/K
#   kco  <- 2.09476e5  # ppm, US Standard Atmosphere
#   
#   vc <- kc25*exp(dhac*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
#   vo <- ko25*exp(dhao*(temp - 25.0)/(298.15*kR*(temp + 273.15)))
#   k  <- vc*(1 + kco*(1e-6)*patm/vo)
#   
#   return(k)
#   
# }
# 
# # $3. calculate Gstar (CO2 compensation point) in Pa
# cal_gstar_gepisat <- function( temp ) {
#   #-----------------------------------------------------------------------
#   # Input:    float, air temperature, degrees C (tc)
#   # Output:   float, gamma-star, Pa (gs)
#   # Features: Returns the temperature-dependent photorespiratory 
#   #           compensation point, Gamma star (Pascals), based on constants 
#   #           derived from Bernacchi et al. (2001) study.
#   # Ref:      Bernacchi et al. (2001), Improved temperature response 
#   #           functions for models of Rubisco-limited photosynthesis, 
#   #           Plant, Cell and Environment, 24, 253--259.
#   #-----------------------------------------------------------------------
#   
#   # Define constants
#   gs25 <- 4.220    # Pa, assuming 25 deg C & 98.716 kPa)
#   dha  <- 37830    # J/mol
#   kR   <- 8.3145   # J/mol/K
#   
#   gs <- gs25 * exp( dha * ( temp - 25.0 ) / ( 298.15 * kR * ( temp + 273.15 ) ) )
#   
#   return( gs )
#   
# }
# 
# # $4. conver CO2 from ppm to Pa
# cal_co2_to_ca <- function( co2, patm ){
#   #-----------------------------------------------------------------------
#   # Input:    - float, annual atm. CO2, ppm (co2)
#   #           - float, monthly atm. pressure, Pa (patm)
#   # Output:   - ca in units of Pa
#   # Features: Converts ca (ambient CO2) from ppm to Pa.
#   #-----------------------------------------------------------------------
#   
#   ca   <- ( 1.e-6 ) * co2 * patm         # Pa, atms. CO2
#   return( ca )
#   
# }
# 
# # $5. calculate the fraction of absorbed PAR
# cal_fapar <- function(LAI){
#   #-----------------------------------------------------------------------
#   ## Input:
#   # LAI: leaf area index, dimensionless
#   
#   ## Output:
#   # fAPAR: the fraction of the absorbed PAR to the incident PAR, dimensionless
#   
#   ## Features: calculate the fraction of absorbed PAR based on given LAI.
#   #-----------------------------------------------------------------------
#   
#   fapar <- 1-exp(-0.5*LAI)
#   
#   return(fapar)
#   
# }
# 
# # $6. calculate vapour pressure deficit
# cal_vpd <- function(RH,Ta){
#   #-----------------------------------------------------------------------
#   ## Input:
#   # RH: relative humidity, percentage
#   # Ta: air temperature, degree C
#   
#   ## Output:
#   # VPD: vapour pressure deficit, kPa
#   
#   ## Features: calculate the vapour pressure defict based on given air temperature and relative humidity.
#   #-----------------------------------------------------------------------
#   
#   VPD  <- 0.611*exp(17.502*Ta/(Ta+240.97))*(1-RH/100)
#   return(VPD)
#   
# }
# 
# # $7. calculate the light use efficiency based on given environmental factors
# cal_lue<-function(Ta,VPD,CO2,elv=NA,patm=NA){
#   
#   #-----------------------------------------------------------------------
#   ## Input:
#   # Ta: air temperature, degree C
#   # PPFD: photosynthetic phtotn flux density, mol photon/m2
#   # VPD: vapour pressure deficit, kPa
#   # CO2: CO2 concentration, ppm
#   # elv: elevation, m
#   # patm: air pressure, kPa
#   
#   ## Output:
#   # LUE: the light use efficiency
#   
#   ## Features: calculate the light use efficiency based on given environmental factors using Pmodel
#   #-----------------------------------------------------------------------
#   
#   Tc.deg_C<-Ta 
#   VPD.kPa  <-VPD
#   CO2.ppm <- CO2 
#   elv <- elv
#   if (identical(NA, elv) && identical(NA, patm)) {
#     rlang::abort("Aborted. Provide either elevation (arugment elv) or atmospheric pressure (argument patm).")
#   }
#   else if (!identical(NA, elv) && identical(NA, patm)) {
#     rlang::warn("Atmospheric pressure (patm) not provided. Calculating it as a function of elevation (elv), assuming standard atmosphere (101325 Pa at sea level).")
#     patm <- calc_patm(elv)
#   } else {
#     patm <- patm*1000 # convert kPa to Pa
#   }
#   
#   # define constant
#   beta <- 146# the ratio of cost factor b to a at reference temperature
#   c <- 0.41# the cost factor of maintaining Jmax
#   
#   # instrinsic quantum yield, based on Cozettle et al.1998, unit: g C/ mol photon
#   # Bernacchi et al. (2001)
#   phi <- (0.352+0.021*Tc.deg_C-3.4*10^(-4)*(Tc.deg_C)^2)/8 # the tempereture-dependence of instrinsic quantum yield of C3
#   maxQE <- phi*12 # convert unit from mol C/ mol photon to g C/ mol photon
#   
#   # light use efficiency
#   K <- cal_k(Tc.deg_C,patm) # the effective Michaelis-Menten coefficient of Rubisco, Pa
#   Gstar <- cal_gstar_gepisat(Tc.deg_C) # photorespiratory compensation point, Pa
#   f1 <- exp(-0.0227*(Tc.deg_C-25)) # the viscosity of water relative to its value at 25藲C
#   ca <- cal_co2_to_ca(CO2.ppm,patm) # ambient CO2 partical pressure, Pa
#   
#   m <- (ca - Gstar)/(ca + 2*Gstar + 3*Gstar*sqrt(1.6*VPD.kPa*1000*f1/(K + Gstar)/(beta)))
#   M <- m*sqrt(1-(c/m)^(2/3))
#   
#   LUE <- M*maxQE # LUE controled by Vcmax, Jmax and instrinsic quantum yield
#   
#   return(LUE)
#   
# }
# 
# # $8. calculate gross primary productivity based on given environmental factors 
# cal_gpp<-function(fAPAR,Ta,PPFD,VPD,CO2,elv=NA,patm=NA){
#   #-----------------------------------------------------------------------
#   ## Input:
#   # fAPAR: the fraction of absorbed PAR to the incident PAR, dimensionless
#   # Ta: air temperature, degree C
#   # PPFD: photosynthetic phtotn flux density, mol photon/m2
#   # VPD: vapour pressure deficit, kPa
#   # CO2: CO2 concentration, ppm
#   # elv: elevation, m
#   # patm: air pressure, kPa
#   
#   ## Output:
#   # GPP: gross primary productivity, g C/m2
#   
#   ## Features: calculate gross primary productivity
#   #-----------------------------------------------------------------------
#   
#   LUE <- cal_lue(Ta = Ta,VPD = VPD,CO2 = CO2,elv = elv,patm = patm)
#   Iabs <- fAPAR*PPFD
#   GPP <- LUE*Iabs
#   
#   return(GPP)
# }
# 
# # $9. calculate aboveground biomass
# cal_ab<-function(GPP_ac,alpha){
#   #-----------------------------------------------------------------------
#   ## Input:
#   # GPP_ac: accumulated gross primary productivity over growing season, g C/m2
#   # alpha: moisture index, the ratio of AET to EET, dimensionless
#   
#   ## Output:
#   # AB: aboveground biomass when harvest, g mass/ m2
#   
#   ## Features: calculate the aboveground biomass based on given alpha and GPP_ac
#   #-----------------------------------------------------------------------
#   
#   # define constant
#   f_tb <- 0.5 # the ration of total biomass to GPP
#   f_C_to_mass <- 2.5 # the conversion coefficient from g C/m2 to g dry mass/m2
#   f_alpha_h <- -0.11 # the sensitivity coefficient of root ratio (root biomass/ total biomass) to alpha
#   
#   # the root ratio when harvest
#   f_root_h <- f_alpha_h*alpha+0.29
#   # the total biomass, g dry mass/m2
#   TB <- GPP_ac*f_tb*f_C_to_mass
#   # aboveground biomass, g dry mass/m2
#   AB<-TB*(1-f_root_h)
#   
#   return(AB)
#   
# }
# 
# # calculate yield based on given aboveground biomass and the amount of nitrogen
# cal_yield <- function(AB){
#   #-----------------------------------------------------------------------
#   ## Input:
#   # AB: aboveground biomass when harvest, g dry mass/ m2
#   
#   ## Output:
#   # yield: yield, g dry mass/m2
#   
#   ## Features: based on given aboveground biomass and the amount of nitrogen
#   #-----------------------------------------------------------------------
#   
#   yield <- 1138.4*(1-exp(-0.00084*AB))-67.5
#   
#   return(yield)
# }
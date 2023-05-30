## Original code Giulia Mengoli and Victor Flo

#' Invokes a P-model function call for sub-daily estimations accounting for acclimation
#'
#' R implementation of the P-model and its 
#' corollary predictions (Prentice et al., 2014; Han et al., 2017; Mengoli et al. 2022).
#'
#' @param TIMESTAMP timestamp in a ymd_hms format (as obtained from the package lubridate)
#' @param tc Temperature, relevant for photosynthesis (deg C). Numeric vector.
#' @param vpd Vapour pressure deficit (Pa)
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param fapar (Optional) Fraction of absorbed photosynthetically active
#'  radiation (unitless, defaults to \code{NA})
#' @param LAI (Optional) Leaf Area Index (m2 leaf m-2 ground, defaults to \code{NA})
#' @param ppfd Incident photosynthetic photon flux density 
#'  (mol m-2 d-1, defaults to \code{NA}). Note that the units of 
#'  \code{ppfd} (per area and per time) determine the units of outputs 
#'  \code{lue}, \code{gpp}, \code{vcmax}, and \code{rd}. For example, 
#'  if \code{ppfd} is provided in units of mol m-2 month-1, then
#'  respective output variables are returned as per unit months.
#' @param u (Optional) Wind speed (m s-1, defaults to \code{NA})
#' @param ustar (Optional) Wind friction velocity (m s-1, defaults to \code{NA})
#' @param sw_in (optional) needed when upscaling_method is "max_rad". 
#' Numeric vector that should be provided in W m-2.
#' @param patm Atmospheric pressure (Pa). When provided, overrides
#'  \code{elv}, otherwise \code{patm} is calculated using standard
#'  atmosphere (101325 Pa), corrected for elevation (argument \code{elv}),
#'  using the function \link{calc_patm}.
#' @param elv Elevation above sea-level (m.a.s.l.). Is used only for 
#'  calculating atmospheric pressure (using standard atmosphere (101325 Pa),
#'  corrected for elevation (argument \code{elv}), using the function
#' \link{calc_patm}), if argument \code{patm} is not provided. If argument
#' \code{patm} is provided, \code{elv} is overridden.
#' @param z (Optional) Wind speed measurement height (m, defaults to \code{NA}).
#' @param d (Optional) numeric, leaf characteristic width (m, defaults to \code{NA})
#' @param netrad (optional) A numeric vector containing Net radiation in W m-2.
#' @param kphio Apparent quantum yield efficiency (unitless). Defaults to
#'  0.081785 for \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, 
#'  do_soilmstress=FALSE}, 0.087182 for \code{method_jmaxlim="wang17",
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}, and 0.049977 for 
#'  \code{method_jmaxlim="wang17", do_ftemp_kphio=FALSE, do_soilmstress=FALSE},
#'  corresponding to the empirically fitted value as presented in Stocker et al.
#'  (2019) Geosci. Model Dev. for model setup 'BRC', 'FULL', and 'ORG' 
#'  respectively, corresponding to \eqn{(a_L b_L)/4} in 
#'  Eq.20 in Stocker et al. (2020) for C3 photosynthesis. For C4 photosynthesis
#'  (\code{c4 = TRUE}), \code{kphio} defaults to 1.0, corresponding to the 
#'   parametrisation by  Cai & Prentice (2020).
#' @param beta Unit cost ratio. Defaults to 146.0 (see Stocker et al., 2019) for
#'   C3 plants and 146/9 for C4 plants.
#' @param c_cost numeric if c_cost is manually added.
#' @param soilm (Optional, used only if \code{do_soilmstress==TRUE}) Relative 
#'  soil moisture as a fraction of field capacity (unitless). Defaults to 1.0 
#'  (no soil moisture stress). This information is used to calculate
#'  an empirical soil moisture stress factor (\link{calc_soilmstress}) whereby
#'  the sensitivity is determined by average aridity, defined by the local 
#'  annual mean ratio of actual over potential evapotranspiration, supplied by
#'  argument \code{meanalpha}.
#' @param meanalpha (Optional, used only if \code{do_soilmstress==TRUE}) Local 
#'  annual mean ratio of actual over potential evapotranspiration, measure for 
#'  average aridity. Defaults to 1.0. Only scalar numbers are accepted. If 
#'  a vector is provided, only the first element will be used.
#' @param apar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) 
#'  Parameter determining the sensitivity of the empirical soil moisture stress 
#'  function. Defaults to 0.0, the empirically fitted value as presented in 
#'  Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' 
#'  (corresponding to a setup with \code{method_jmaxlim="wang17", 
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param bpar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) 
#'  Parameter determining the sensitivity of the empirical soil moisture stress
#'  function. Defaults to 0.7330, the empirically fitted value as presented in 
#'  Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' 
#'  (corresponding to a setup with \code{method_jmaxlim="wang17", 
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param c4 (Optional) A logical value specifying whether the C3 or C4
#'  photosynthetic pathway is followed.Defaults to \code{FALSE}. If \code{TRUE},
#'  the leaf-internal CO2 concentration is still estimated using beta but
#'  \eqn{m} (returned variable \code{mj}) tends to 1, and \eqn{m'} tends to
#'  0.669 (with \code{c = 0.41}) to represent CO2 concentrations within the leaf.
#'  With \code{do_ftemp_kphio = TRUE}, a C4-specific temperature dependence of
#'  the quantum yield efficiency is used (see \link{ftemp_kphio}).
#' @param method_jmaxlim (Optional) A character string specifying which method 
#'  is to be used for factoring in Jmax limitation. Defaults to \code{"wang17"},
#'  based on Wang Han et al. 2017 Nature Plants and (Smith 1937). Available is 
#'  also \code{"smith19"}, following the method by Smith et al., 2019 Ecology 
#'  Letters, and \code{"none"} for ignoring effects of Jmax limitation.
#' @param do_ftemp_kphio (Optional) A logical specifying whether 
#'  temperature-dependence of quantum yield efficiency is used. See \link{ftemp_kphio}
#'  for details. Defaults to \code{TRUE}. Only scalar numbers are accepted. If 
#'  a vector is provided, only the first element will be used.
#' @param do_soilmstress (Optional) A logical specifying whether an empirical 
#' soil moisture stress factor is to be applied to down-scale light use 
#' efficiency (and only light use efficiency). Defaults to \code{FALSE}.
#' @param returnvar (Optional) A character string of vector of character strings
#'  specifying which variables are to be returned (see return below).
#' @param verbose Logical, defines whether verbose messages are printed. 
#'  Defaults to \code{FALSE}.
#' @param upscaling_method String defining the method applied in the sub-daily down scaling. 
#'  Accepted values are "noon", "daily" and "max_rad". "noon" approach computes the 
#'  down scaling using the average of conditions around midday (hour_reference_t +/- nr_window).
#'  "daily" approach uses average daytime conditions. "max_rad" uses average of daytime 
#'  conditions around the point of maximum radiation.
#' @param hour_reference_T numeric from 0 to 23. Reference time for the upscaling process.
#' @param gap_method character string of the method used to do gapfilling. Accepted values are "linear" or "continous".
#' @param acclim_days numeric with the number of days to use for acclimation.
#'  
#'
#' @return A named list of numeric values (including temperature and pressure 
#' dependent parameters of the photosynthesis model, P-model predictions, 
#' including all its corollary). This includes :
#' 
#' \itemize{
#'  \item \code{ca}: Ambient CO2 expressed as partial pressure (Pa)
#'  
#'  \item \code{gammastar}: Photorespiratory compensation point \eqn{\Gamma*}, 
#'   (Pa), see \link{calc_gammastar}.
#'  
#'  \item \code{kmm}: Michaelis-Menten coefficient \eqn{K} for photosynthesis 
#'  (Pa), see \link{calc_kmm}.
#'  
#'  \item \code{ns_star}: Change in the viscosity of water, relative to its 
#'   value at 25 deg C (unitless).
#'   \deqn{\eta* = \eta(T) / \eta(25 deg C)}
#'   This is used to scale the unit cost of transpiration. 
#'   Calculated following Huber et al. (2009).
#'  
#'  \item \code{chi}: Optimal ratio of leaf internal to ambient CO2 (unitless). 
#'   Derived following Prentice et al.(2014) as:
#'  \deqn{
#'   \chi = \Gamma* / ca + (1- \Gamma* / ca) \xi / (\xi + \sqrt D )
#'   }
#'   with
#'   \deqn{
#'    \xi = \sqrt (\beta (K+ \Gamma*) / (1.6 \eta*))
#'   }
#'   \eqn{\beta} is given by argument \code{beta}, \eqn{K} is 
#'   \code{kmm} (see \link{calc_kmm}), \eqn{\Gamma*} is 
#'   \code{gammastar} (see \link{calc_gammastar}). \eqn{\eta*} is \code{ns_star}.
#'   \eqn{D} is the vapour pressure deficit (argument \code{vpd}), \eqn{ca} is 
#'   the ambient CO2 partial pressure in Pa (\code{ca}).
#'   
#'   \item \code{ci}: Leaf-internal CO2 partial pressure (Pa), calculated as \eqn{(\chi ca)}.
#'   
#'   \item \code{lue}: Light use efficiency (g C / mol photons), calculated as
#'                         \deqn{
#'                              LUE = \phi(T) \phi0 m' Mc
#'                         }
#'                         where \eqn{\phi(T)} is the temperature-dependent quantum yield efficiency modifier
#'                         (\link{ftemp_kphio}) if \code{do_ftemp_kphio==TRUE}, and 1 otherwise. \eqn{\phi 0}
#'                         is given by argument \code{kphio}.
#'                         \eqn{m'=m} if \code{method_jmaxlim=="none"}, otherwise
#'                         \deqn{
#'                                m' = m \sqrt( 1 - (c/m)^(2/3) )
#'                         }
#'                         with \eqn{c=0.41} (Wang et al., 2017) if \code{method_jmaxlim=="wang17"}. \eqn{Mc} is
#'                         the molecular mass of C (12.0107 g mol-1). \eqn{m} is given returned variable \code{mj}.
#'                         If \code{do_soilmstress==TRUE}, \eqn{LUE} is multiplied with a soil moisture stress factor,
#'                         calculated with \link{calc_soilmstress}.
#'         \item \code{mj}: Factor in the light-limited assimilation rate function, given by
#'                         \deqn{
#'                             m = (ci - \Gamma*) / (ci + 2 \Gamma*)
#'                        }
#'                        where \eqn{\Gamma*} is given by \code{calc_gammastar}.
#'         \item \code{mc}: Factor in the Rubisco-limited assimilation rate function, given by
#'                         \deqn{
#'                             mc = (ci - \Gamma*) / (ci + K)
#'                        }
#'                        where \eqn{K} is given by \code{calc_kmm}.
#'         \item \code{gpp}: Gross primary production (g C m-2), calculated as
#'                        \deqn{
#'                            GPP = Iabs LUE
#'                        }
#'                        where \eqn{Iabs} is given by \code{fapar*ppfd} (arguments), and is
#'                        \code{NA} if \code{fapar==NA} or \code{ppfd==NA}. Note that \code{gpp} scales with
#'                        absorbed light. Thus, its units depend on the units in which \code{ppfd} is given.
#'         \item \code{iwue}: Intrinsic water use efficiency (iWUE, Pa), calculated as
#'                        \deqn{
#'                              iWUE = ca (1-\chi)/(1.6)
#'                        }
#'         \item \code{gs}: Stomatal conductance (gs, in mol C m-2 Pa-1), calculated as
#'                        \deqn{
#'                             gs = A / (ca (1-\chi))
#'                        }
#'                        where \eqn{A} is \code{gpp}\eqn{/Mc}.
#'         \item \code{vcmax}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) at growth temperature (argument
#'                       \code{tc}), calculated as
#'                       \deqn{
#'                            Vcmax = \phi(T) \phi0 Iabs n
#'                       }
#'                       where \eqn{n} is given by \eqn{n=m'/mc}.
#'         \item \code{vcmax25}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) normalised to 25 deg C
#'                      following a modified Arrhenius equation, calculated as \eqn{Vcmax25 = Vcmax / fv},
#'                      where \eqn{fv} is the instantaneous temperature response by Vcmax and is implemented
#'                      by function \link{ftemp_inst_vcmax}.
#'         \item \code{jmax}: The maximum rate of RuBP regeneration () at growth temperature (argument
#'                       \code{tc}), calculated using
#'                       \deqn{
#'                            A_J = A_C
#'                       }
#'         \item \code{rd}: Dark respiration \eqn{Rd} (mol C m-2), calculated as
#'                      \deqn{
#'                          Rd = b0 Vcmax (fr / fv)
#'                      }
#'                      where \eqn{b0} is a constant and set to 0.015 (Atkin et al., 2015), \eqn{fv} is the
#'                      instantaneous temperature response by Vcmax and is implemented by function
#'                      \link{ftemp_inst_vcmax}, and \eqn{fr} is the instantaneous temperature response
#'                      of dark respiration following Heskel et al. (2016) and is implemented by function
#'                      \link{ftemp_inst_rd}.
#' }
#'
#' Additional variables are contained in the returned list if argument \code{method_jmaxlim=="smith19"}
#' \itemize{
#'  \item \code{omega}: Term corresponding to \eqn{\omega}, defined by Eq. 16 in
#'   Smith et al. (2019), and Eq. E19 in Stocker et al. (2019).
#'   
#'  \item \code{omega_star}: Term corresponding to \eqn{\omega^\ast}, defined by
#'   Eq. 18 in Smith et al. (2019), and Eq. E21 in Stocker et al. (2019).
#'  }patm
#'
#' @references  
#'  Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature response func-tions  of  parameters
#'  required  to  model  RuBP-limited  photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#'
#   Cai, W., and Prentice, I. C.: Recent trends in gross primary production 
#'  and their drivers: analysis and modelling at flux-site and global scales,
#'  Environ. Res. Lett. 15 124050 https://doi.org/10.1088/1748-9326/abc64e, 2020
#
#'  Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J.,
#'  Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K.,
#'  Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response
#'  of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#'  113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#'
#'  Huber,  M.  L.,  Perkins,  R.  A.,  Laesecke,  A.,  Friend,  D.  G.,  Sengers,  J.  V.,  Assael,M. J.,
#'  Metaxa, I. N., Vogel, E., Mares, R., and Miyagawa, K.:  New international formulation for the viscosity
#'  of H2O, Journal of Physical and Chemical ReferenceData, 38, 101–125, 2009
#'
#'  Mengoli, G., Agustí-Panareda, A., Boussetta, S., Harrison, S. P., Trotta, C., 
#'  & Prentice, I. C. (2022). Ecosystem photosynthesis in land-surface models: A first-principles approach 
#'  incorporating acclimation. Journal of Advances in Modeling Earth Systems, 14,
#'   e2021MS002767. https://doi.org/10.1029/2021MS002767 
#'  
#'  Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancing the costs
#'  of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology,
#'  Ecology  Letters,  17,  82–91,  10.1111/ele.12211,http://dx.doi.org/10.1111/ele.12211, 2014.
#'
#'  Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J.,
#'  and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.
#'  Atkin, O. K., et al.:  Global variability in leaf respiration in relation to climate, plant func-tional
#'  types and leaf traits, New Phytologist, 206, 614–636, doi:10.1111/nph.13253,
#'  https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.13253.
#'
#'  Smith, N. G., Keenan, T. F., Colin Prentice, I. , Wang, H. , Wright, I. J., Niinemets, U. , Crous, K. Y.,
#'  Domingues, T. F., Guerrieri, R. , Yoko Ishida, F. , Kattge, J. , Kruger, E. L., Maire, V. , Rogers, A. ,
#'  Serbin, S. P., Tarvainen, L. , Togashi, H. F., Townsend, P. A., Wang, M. , Weerasinghe, L. K. and Zhou, S.
#'  (2019), Global photosynthetic capacity is optimized to the environment. Ecol Lett, 22: 506-517.
#'  doi:10.1111/ele.13210
#'
#'  Stocker, B. et al. Geoscientific Model Development Discussions (in prep.)
#'
#' @export
#'
#' @examples \dontrun{
#'  rpmodel(
#'   tc = 20,
#'   vpd = 1000,
#'   co2 = 400,
#'   ppfd = 30,
#'   elv = 0)
#' }
#'
rpmodel_subdaily <- function(
    TIMESTAMP, tc, vpd, co2, fapar = NA, LAI = NA, ppfd, u = NA, ustar = NA,#wind speed in m s^-1
    canopy_height=NA, sw_in = NA, patm = NA, elv = NA, z = NA, leafwidth = NA, netrad = NA,
    kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.087182, 0.081785), 0.049977),
    beta = 146.0, c_cost = 0.41, soilm = 1.0, meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.73300,
    c4 = FALSE, method_jmaxlim = "wang17",do_ftemp_kphio = TRUE, do_soilmstress = FALSE,
    do_leaftemp = FALSE, gb_method = "Su_2001", do_acclimation = FALSE, epsleaf = 0.96, #thermal absorptivity of the leaf
    energy_params = list(
      ste_bolz = 5.67e-8, #W m^-2 K^-4
      cpm = 75.38, #J mol^-1 ºC-1
      J_to_mol = 4.6, #Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
      frac_PAR = 0.5, #Fraction of incoming solar irradiance that is photosynthetically active radiation (PAR
      fanir = 0.35 #Fraction of NIR absorbed
    ), returnvar = NULL, verbose = FALSE,
    upscaling_method = c("noon","daily","max_rad"), hour_reference_T = 12, gap_method = "linear",
    acclim_days = 15, weighted_accl = TRUE
){
  
  #---- Fixed parameters--------------------------------------------------------
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  #'
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  
  #parameters
  # epsleaf = energy_params["epsleaf"] %>% as.numeric
  sigma = energy_params["ste_bolz"]%>% as.numeric
  cpm = energy_params["cpm"]%>% as.numeric #molar heat capacity of water 
  # kfFEC = energy_params["kfFEC"]%>% as.numeric
  fanir = energy_params["fanir"]%>% as.numeric
  J_to_mol = energy_params["J_to_mol"]%>% as.numeric
  frac_PAR = energy_params["frac_PAR"]%>% as.numeric
  lat_heat = 2230 #J g^-1
  mol_gas_const = 8.3144621 #J mol^-1 K^-1
  mol_mas_wv = 18.01528 #g mol-1
  mol_mas_dry_air = 28.9652 #g mol-1
  CP <- 1005 # specific heat capacity of air at constant pressure (J kg-1 K-1)
  spe_gas_const = mol_gas_const/mol_mas_wv #J g^-1 K^-1
  
  #---- soil moisture stress as a function of soil moisture and mean alpha -----
  if (do_soilmstress) {
    if (length(AI) > 1){
      warning("Argument 'AI' has length > 1. Only the first element is used.")
      AI <- AI[1]
    }
    soilmstress <- calc_soilmstress(soilm, AI)
  }
  else {
    soilmstress <- 1.0
  }
  
  #---- check for negative ppfd ----
  if(any(ppfd<0)){
    warning("Some values of PPFD are negative. Those PPFD < 0 will be set to 0.")
    ppfd[ppfd<0] <- 0
  }
  
  # 1.0 Calculate P model without acclimation
  # tibble(TIMESTAMP,tc, vpd, co2, fapar, LAI, ppfd, u, ustar, canopy_height, sw_in, netrad, patm, meanalpha) %>% 
  #   split(seq(nrow(.))) %>%
  #   purrr::map_df(function(x){
    #   res <- rpmodel(x$tc, x$vpd, x$co2, x$fapar, x$LAI,
    #                 x$ppfd, x$u, x$ustar, x$canopy_height,x$sw_in, x$patm, 
    #                 unique(elv), unique(z),x$netrad, kphio, beta, c_cost,
    #                 soilm, x$meanalpha, apar_soilm, bpar_soilm, c4,
    #                 method_jmaxlim, do_ftemp_kphio, do_soilmstress, do_leaftemp = FALSE,
    #                 verbose = verbose, energy_params = energy_params)%>% 
    #   as_tibble()
    # return(res)
    # }) 
  
  rpmodel(tc, vpd, co2, fapar, LAI,
          ppfd, u, ustar, canopy_height, sw_in, patm, 
          elv, z, leafwidth, netrad, kphio, beta, c_cost,
          soilm, AI, c4, method_jmaxlim, do_ftemp_kphio, 
          do_soilmstress, do_leaftemp = FALSE,
          gb_method = gb_method, epsleaf = epsleaf,
          energy_params = energy_params, verbose = verbose)%>% 
    as_tibble()-> df_Or
  
  #DO ACCLIMATION?
  if(!do_acclimation){
    # print("No acclimation was calculated")
    return(df_Or)
  }else{ 
  
  
    # 2.0 Create df for upscaling
    dfIn <- tibble(TIMESTAMP = TIMESTAMP,
                   YEAR = lubridate::year(TIMESTAMP),
                   MONTH = lubridate::month(TIMESTAMP),
                   DAY = lubridate::day(TIMESTAMP),
                   HOUR = lubridate::hour(TIMESTAMP),
                   MINUTE = lubridate::minute(TIMESTAMP),
                   fapar = fapar,
                   tc = tc,
                   ppfd = ppfd,
                   vpd = vpd,
                   co2 = co2,
                   u = u,
                   ustar = ustar,
                   z = z,
                   netrad = netrad,
                   LAI = LAI,
                   canopy_height = canopy_height,
                   patm = patm,
                   epsleaf = epsleaf,
                   soilm = soilm,
                   soilmstress = soilmstress)
    
    # 2.1 apply the dailyUpscaling function 
    dataDaily <- dailyUpscaling(df = dfIn, 
                                  nrWindow = 1, 
                                  hour_reference_T = hour_reference_T, 
                                  upscaling_method = upscaling_method)
    
    # 3.0 apply running mean
    # dataDaily = runningMean(data_x_running_mean_t = dataDaily, daily_window = acclim_days)
    if(!weighted_accl){
    dataDaily = lapply(dataDaily[,-c(1:6)],function(x)runner::mean_run(x= x,k=acclim_days)) %>% as.data.frame() %>% cbind(dataDaily[,c(1:6)])
    }else{
    dataDaily = lapply(dataDaily[,-c(1:6)] %>% dplyr::select_if(function(x) any(!is.na(x))) ,function(y){
      # berryFunctions::movAv(dat=y,width=acclim_days,weights = rep(1, acclim_days)- seq(0,1, length.out = acclim_days))})%>% 
      # as.data.frame() %>% cbind(dataDaily[,c(1:6)])
      # TTR::EMA(x=y,n=acclim_days)})%>% 
      # as.data.frame() %>% cbind(dataDaily[,c(1:6)])
      pracma::movavg(x=y,n=acclim_days,type="e")})%>% 
      as.data.frame() %>% cbind(dataDaily[,c(1:6)]) %>%cbind(dataDaily %>% dplyr::select_if(function(x) all(is.na(x))))
    }
    
    # 4.0 Calculate P-model on daily upscaling values
    tibble(dataDaily) %>% 
      split(seq(nrow(.))) %>%
      purrr::map_df(function(x){
        # print(paste("acclimation of", x$TIMESTAMP))
        res <- rpmodel(x$tc, x$vpd, x$co2, x$fapar, x$LAI,
                       x$ppfd, x$u, x$ustar, x$canopy_height,x$sw_in, x$patm, 
                       unique(elv), unique(z), leafwidth, x$netrad, kphio, beta, c_cost,
                       x$soilm, AI, c4,
                       method_jmaxlim, do_ftemp_kphio, do_soilmstress, 
                       do_leaftemp = do_leaftemp, gb_method = gb_method, 
                       verbose = verbose, epsleaf=x$epsleaf, energy_params = energy_params)%>% 
          as_tibble()
        return(res)
      }) -> df_Mm
    
    
    # 5.0 Downscale from daily to subdaily
    dataDaily = dataDaily %>% cbind(df_Mm) 
    
    DF <- dfIn  %>% 
      cbind(df_Or) %>% 
      left_join(dataDaily %>%
                  dplyr::select(where(~sum(!is.na(.x)) > 0))%>% 
                  dplyr::select(-c("YEAR","MONTH","DAY","HOUR","MINUTE")) %>% 
                  rename_with( ~ paste0(.x, "_opt"), where(is.numeric)),
                by = "TIMESTAMP")
    
    
    # 5.1 GAP- FILLING
    DF <- lapply(DF, function(y) {
      if (is.numeric(y)) { 
        approx(as.numeric(DF$TIMESTAMP), y, method = gap_method, 
               xout = as.numeric(DF$TIMESTAMP), rule = 2)$y
        }else{
          y
        }
      }
      ) %>% bind_cols()
    
  
    
    
    # 6.1 Optimal leaf temperature
    if(do_leaftemp){
    tcleaf_new <- DF %>% 
      split(seq(nrow(.))) %>%
      purrr::map_df(function(x){
        tryCatch({
          # print(x$TIMESTAMP)
          rho <- calc_air_density(x$patm, x$tc)  # air density (kg m-3)
          es = exp(34.494-4924.99/(x$tc+237.1))/((x$tc+105)^1.57)
          ea = es - x$vpd
          patm = x$patm
          tk = x$tc+273.15
          tcleaf_dew <- tk/(1-17.27^(-1)*log(ea/es))-273.15
          # Td = (243.5 * log(ea/es) + (17.67 * x$tc)/(x$tc + 243.5)) / 
          #   (17.67 - log(ea/es) - (17.67 * x$tc)/(x$tc + 243.5))
          # tcleaf = x$tc
          # "2015-05-09 01:00:00 UTC"
          #################################################
          tcleaf_new <- optimise(function(tcleaf_root){
            # print(tcleaf_root)
            tkleaf = tcleaf_root+273.15
            # ei = es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tkleaf)) #assuming saturation within the leaf (We don't apply psi_leaf correction es(T)exp(ψleaf V/(RT)))
            ei = exp(34.494-4924.99/(tcleaf_root+237.1))/((tcleaf_root+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
            vpd_new = (ei - ea)
            vpd_new = ifelse(vpd_new<0,0,vpd_new)

            EB <- energy_balance(tcleaf=tcleaf_root, tcleaf_opt = x$tcleaf_opt, 
                                 vpd = vpd_new, ppfd = x$ppfd, ppfd_opt = x$ppfd_opt, 
                                 fapar = x$fapar, fapar_opt = x$fapar_opt, ca = x$ca, 
                                 ca_opt = x$ca_opt, xi = x$xi, xiPa = x$xiPa, 
                                 patm = x$patm, ns_star_opt = x$ns_star_opt, 
                                 gammastar = x$gammastar, gammastar_opt = x$gammastar_opt, 
                                 kmm = x$kmm, kmm_opt = x$kmm_opt, kphio = kphio,
                                 soilmstress = soilmstress, method_jmaxlim = method_jmaxlim,
                                 c4 = c4, rd_to_vcmax = rd_to_vcmax, 
                                 beta = beta, c_cost = c_cost, 
                                 u = x$u, canopy_height = x$canopy_height,
                                 tc = x$tc, tk = tk, tkleaf = tkleaf, 
                                 z = x$z, LAI = x$LAI, ustar = x$ustar, netrad = x$netrad,
                                 mol_gas_const =  mol_gas_const, J_to_mol = J_to_mol, 
                                 lat_heat = lat_heat, mol_mas_wv = mol_mas_wv, 
                                 sigma = sigma, cpm = cpm, CP = CP, rho = rho,
                                 epsleaf = x$epsleaf, epssky = epssky, frac_PAR = frac_PAR, 
                                 fanir = fanir, ei = ei, ea = ea, gb_method = gb_method, 
                                 leafwidth = leafwidth)
            if(is.na(x$netrad)){
              Rnet = EB$Rnet
            }else{
              Rnet = x$netrad
            }
            #final balance
            (Rnet - EB$Qc - EB$lE)^2
          },
          interval=c(tcleaf_dew, x$tc+30))$minimum
          ##################################################
          
          tkleaf = tcleaf_new+273.15
          # ei = es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tkleaf)) #assuming saturation within the leaf (We don't apply psi_leaf correction es(T)exp(ψleaf V/(RT)))
          ei = exp(34.494-4924.99/(tcleaf_new+237.1))/((tcleaf_new+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
          vpd_new = (ei - ea)
          vpd_new = ifelse(vpd_new<0,0,vpd_new)
          EB <- energy_balance(tcleaf=tcleaf_new, tcleaf_opt = x$tcleaf_opt, 
                               vpd = vpd_new, ppfd = x$ppfd, ppfd_opt = x$ppfd_opt, 
                               fapar = x$fapar, fapar_opt = x$fapar_opt, ca = x$ca, 
                               ca_opt = x$ca_opt, xi = x$xi, xiPa = x$xiPa, 
                               patm = x$patm, ns_star_opt = x$ns_star_opt, 
                               gammastar = x$gammastar, gammastar_opt = x$gammastar_opt, 
                               kmm = x$kmm, kmm_opt = x$kmm_opt, kphio = kphio,
                               soilmstress = soilmstress, method_jmaxlim = method_jmaxlim,
                               c4 = c4, rd_to_vcmax = rd_to_vcmax, 
                               beta = beta, c_cost = c_cost, 
                               u = x$u, canopy_height = x$canopy_height,
                               tc = x$tc, tk = tk, tkleaf = tkleaf, 
                               z = x$z, LAI = x$LAI, ustar = x$ustar, netrad = x$netrad,
                               mol_gas_const =  mol_gas_const, J_to_mol = J_to_mol, 
                               lat_heat = lat_heat, mol_mas_wv = mol_mas_wv, 
                               sigma = sigma, cpm = cpm, CP = CP, rho = rho,
                               epsleaf = x$epsleaf, epssky = epssky, frac_PAR = frac_PAR, 
                               fanir = fanir, ei = ei, ea = ea,  gb_method = gb_method,
                               leafwidth = leafwidth)
          
        return(tibble(tcleaf_new=tcleaf_new ,Q_tcleaf_new = 1, Qc = EB$Qc, 
                      tcleaf_dew = tcleaf_dew, Rnet = EB$Rnet, lE = EB$lE, 
                      Qtir = EB$Qtir, Qtirleaf = EB$Qtirleaf, gb = EB$gb, 
                      ust = EB$ust))
          },
        error= function(e){
          return(tibble(tcleaf_new=x$tc, Q_tcleaf_new = 0, tcleaf_dew = tcleaf_dew,
                        Qc = NA, Rnet = NA, le = NA, Qtir = NA, Qtirleaf = NA, 
                        gb = NA, ust = NA))
          }
        )
    }) %>% bind_rows() 
    tcleaf = tcleaf_new$tcleaf_new %>% as.numeric()
    tkleaf = tcleaf+273.15
    es = exp(34.494-4924.99/(tc+237.1))/((tc+105)^1.57)
    ea = es - vpd
    ei = exp(34.494-4924.99/(tcleaf+237.1))/((tcleaf+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
    vpd_new = (ei - ea)
    vpd_new = ifelse(vpd_new<0,0,vpd_new)
    }else{
      tcleaf = tc
      vpd_new = vpd 
    }
    
    # 7 Calculate acclimated subdaily
    out <- rpmodel_jmax_vcmax(tcleaf=tcleaf, tcleaf_opt = DF$tcleaf_opt, vpd = vpd_new, ppfd = DF$ppfd, ppfd_opt = DF$ppfd_opt,
                              fapar = DF$fapar, fapar_opt = DF$fapar_opt, ca = DF$ca, ca_opt = DF$ca_opt, 
                              xi = DF$xi, xiPa = DF$xiPa, patm = DF$patm, ns_star_opt = DF$ns_star_opt, 
                              gammastar = DF$gammastar, gammastar_opt = DF$gammastar_opt, kmm = DF$kmm, kmm_opt = DF$kmm_opt,
                              kphio = kphio, soilmstress = soilmstress, method_jmaxlim = method_jmaxlim, c4 = c4, rd_to_vcmax = rd_to_vcmax,
                              beta = beta, c_cost = c_cost, leafwidth = leafwidth, LAI = LAI)
    
    
    if(do_leaftemp){
      out$Q_tcleaf <-  tcleaf_new$Q_tcleaf_new
      out$Qc <-  tcleaf_new$Qc
      out$Rnet <-  tcleaf_new$Rnet
      out$gb <-  tcleaf_new$gb
      out$Qtir <-  tcleaf_new$Qtir
      out$Qtirleaf <-  tcleaf_new$Qtirleaf
      out$lE <-  tcleaf_new$lE
      out$ust <-  tcleaf_new$ust
      out$tcleaf_dew <- tcleaf_new$tcleaf_dew
    }
    
    return(out)
  }
}




# df = dfIn
# nrWindow = 1
# hour_reference_T = c(10,12,15)
# upscaling_method = "max_rad"
# cicleday = 1
# ciclemonth = 11
# cicleyear = 2014
# hourReference = 12

dailyUpscaling <- function(df = dfIn, nrWindow = 1, hour_reference_T = c(10,12,15), 
                             upscaling_method = "noon") {
  
  #1.0 Header control
  if(!is.character(upscaling_method)){
    stop(
      cat('The upscaling_method should be a character string either "noon", "daily" or "max_rad".')
    )
  }else  if(!upscaling_method %in% c("noon","daily","max_rad")){
    stop(
      cat('The upscaling_method provided is not a valid value. It should be either "noon", "daily" or "max_rad"')
      )
  }
  
  colMandatory = c('YEAR','MONTH','DAY','HOUR','MINUTE', 'TIMESTAMP')
  # if (upscaling_method == "max_rad")
  #   colMandatory = c(colMandatory,"sw_in")
  
  if (!headerControl_dd(df = df, colMandatory = colMandatory))
    stop(headerControl_dd(df = df, colMandatory = colMandatory, showMsg = TRUE))
  
  #2.0 T calculation
  df <- df %>% as.data.frame()
  dfDayT = df[1,]    
  
  
  for (cicleyear in sort(unique(df$YEAR))) {
    for (ciclemonth in seq(1,12)) {
      for (cicleday in seq(1,31)) {
        
        posDay = which(df$YEAR == cicleyear & df$MONTH == ciclemonth &
                         df$DAY == cicleday)
        if (length(posDay) == 0) next
        
        for (hourReference in hour_reference_T) {
          # upscaling_method noon
          if (upscaling_method == "noon") {
            posReference = which(df$YEAR == cicleyear & df$MONTH == ciclemonth & df$DAY == cicleday &
                                   df$HOUR == hourReference & df$MINUTE == 0)
            if (length(posReference) == 0) next
            windowDay = seq(posReference - nrWindow,posReference + nrWindow)
          }
          # upscaling_method daily
          if (upscaling_method == "daily") {
            posReference = which(df$YEAR == cicleyear & df$MONTH == ciclemonth & df$DAY == cicleday &
                                   df$HOUR == 12 & df$MINUTE == 0)
            refdate = lubridate::ymd_hms(paste(cicleyear,ciclemonth,cicleday,"12:00:00"))
            windowDay = which(df$YEAR == cicleyear & df$MONTH == ciclemonth & df$DAY == cicleday)
          }
          # upscaling_method max_rad
          if (upscaling_method == "max_rad") {
            posReferenceMax = which(df$YEAR == cicleyear & df$MONTH == ciclemonth & df$DAY == cicleday &
                                      is.na(df$ppfd) == 0)
            if ( length(posReferenceMax) > 0 ) {
              
              posMaxppfd = which(df$ppfd[posReferenceMax] == max(df$ppfd[posReferenceMax],na.rm = T)[1])
              
              posReference = posReferenceMax[posMaxppfd[1]]
              windowDay = seq(posReference - nrWindow,posReference + nrWindow)
            } else {
              posReference = which(df$YEAR == cicleyear & df$MONTH == ciclemonth & df$DAY == cicleday &
                                     df$HOUR == 12 & df$MINUTE == 0)
              windowDay = NA
            }
          }
          
          
          dfDay = df[1,]
          for (col in colnames(df)) {
            
            if(col != 'TIMESTAMP'){
              
              if (is.na(sum(windowDay)) ) {
                dfDay[1,col] = NA  
                next
              }
            
              tmp = df[windowDay,col]
              if (!is.numeric(tmp)) {
                dfDay[1,col] = NA
                next
              }
            
              
              posNa = which(is.na(tmp) == TRUE)
            
              if (length(posNa) > 0){
                tmp = tmp[-1*posNa]
                }
            
              if (length(tmp) == 0){
                dfDay[1,col] = NA
              } else {
                dfDay[1,col] = mean(tmp,na.rm=TRUE)
              }
              rm(tmp,posNa)
            }else{next}
          }
          
          if (upscaling_method == "daily"){
            dfDay$TIMESTAMP = refdate
          }else{
            dfDay$TIMESTAMP = df$TIMESTAMP[posReference]
            }
          
          dfDayT = rbind(dfDayT,dfDay)  
          
          rm(dfDay,windowDay,posReference)
          
          if (upscaling_method == "daily") break;
          if (upscaling_method == "max_rad") break;
        }
      }
    }
  }
  dfDayT = dfDayT[-1,]
  
  dfDayT = dfDayT[order(dfDayT$TIMESTAMP),]
  
  return(dfDayT)
}




#' Auxiliar function for dailydounscaling function that check if the mandatory variables 
#' exist in the dataset (missing or redundant columns)
#' @param df dataframe to check the colnames
#' @param colMandatory list of mandatory headers
#' @param showMsg Logical. If TRUE, it shows the function messages. Defaults to \code{FALSE}.
#' @return TRUE or FALSE. 

headerControl_dd <- function(df = dfToCheck, colMandatory = listMandatoryToCheck,showMsg = F) {

  # string of mandatory variables, missing
  errorColMissing = c()

  # string of mandatory variables, duplicate
  errorColMultiple = c()

  for (col in colMandatory){
    ckMandatory = which(colnames(df) == col)
    if (length(ckMandatory) == 0) {
      errorColMissing = c(errorColMissing,col)
      next
    }
    if (length(ckMandatory) > 1) {
      errorColMultiple = c(errorColMultiple,col)
      next
    }
  }
  rm(col)
  if (showMsg) {
    cat(sprintf('missing mandatory variables: %d \n',length(errorColMissing)))
    if (length(errorColMissing) > 0)
      cat(sprintf(' (%s)\n',paste(errorColMissing,collapse = ',')))

    cat(sprintf('multiple mandatory variables: %d \n',length(errorColMultiple)))
    if (length(errorColMultiple) > 0)
      cat(sprintf(' (%s)\n',paste(errorColMultiple,collapse = ',')))
  }
  if ( (length(errorColMissing) + length(errorColMultiple)) > 0) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}



# Function to apply the running mean method on the model inputs 
# param: data_x_running_mean_t, dataframe to use for the runningMean 
# param: daily_window, number of days to compute the running mean 
# param: showMsg (if T shows the function messages)
# return: data.frame with runningMean output
# rdname: runningMean


# runningMean <- function(data_x_running_mean_t = df, daily_window = 10) {
#   
#   unique_hour = sort(unique(data_x_running_mean_t$HOUR))
#   
#   data_running_mean_t = data_x_running_mean_t[1,]
#   
#   for (cy_hour_ref in unique_hour) {
#     pos_hour = which(data_x_running_mean_t$HOUR == cy_hour_ref)
#     
#     data_x_running_mean = data_x_running_mean_t[pos_hour,]
#     rm(pos_hour)
#     # window of time 
#     lunghezza_finestra = daily_window - 1# remotion of one day since the computation starts from the 1st
#     data_running_mean = data_x_running_mean[1,]
#     
#     for ( ciclo_inizio_mm in seq(2,nrow(data_x_running_mean)) ) {
#       # definition of positions to compute the running mean by making a counter in reverse
#       posizioni_rm = (ciclo_inizio_mm - lunghezza_finestra):ciclo_inizio_mm
#       
#       pos_meno1 = which(posizioni_rm < 1)
#       if ( length(pos_meno1) > 0 ) posizioni_rm = posizioni_rm[-1*pos_meno1]
#       rm(pos_meno1)
#       
#       # cycle for computing the mean of the dataset's variables
#       media_1 = data_x_running_mean[1,]
#       for ( ciclo_variabili in colnames(media_1) ) {
#         # dataset creation with values to use
#         data_1 = data_x_running_mean[posizioni_rm,ciclo_variabili]
#         if ( !is.numeric(data_1) ) {
#           media_1[ciclo_variabili] = NA
#           next
#         }
#         # NA remotion
#         pos_na = which(is.na(data_1) == 1 )
#         if ( length(pos_na) > 0 ) data_1 = data_1[-1*pos_na]
#         rm(pos_na)
#         # if there are missing values it puts NA, otherwise it computes the mean
#         if ( length(data_1) == 0 ) {
#           media_1[ciclo_variabili] = NA
#         } else {
#           media_1[ciclo_variabili] = mean(data_1)
#         }
#         rm(data_1)
#       }
#       rm(ciclo_variabili)
#       media_1$TIMESTAMP = data_x_running_mean$TIMESTAMP[ciclo_inizio_mm]
#       data_running_mean = rbind(data_running_mean,media_1)
#       rm(media_1)
#       rm(posizioni_rm)
#     }
#     rm(ciclo_inizio_mm)
#     rm(lunghezza_finestra)
#     rm(data_x_running_mean)
#     data_running_mean_t = rbind(data_running_mean_t,data_running_mean)
#   }
#   
#   data_running_mean_t = data_running_mean_t[-1,]
#   
#   # library(lubridate)
#   # data_running_mean_t$YEAR = year(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
#   # data_running_mean_t$MONTH = month(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
#   # data_running_mean_t$DAY = day(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
#   # data_running_mean_t$HOUR = hour(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
#   # data_running_mean_t$MINUTE = minute(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
#   
#   data_running_mean_t = data_running_mean_t[order(data_running_mean_t$TIMESTAMP),]
#   
#   return(data_running_mean_t)
# }



rpmodel_jmax_vcmax <- function(tcleaf, tcleaf_opt, vpd, ppfd, ppfd_opt, fapar, fapar_opt, ca, ca_opt, xi, xiPa, patm, ns_star_opt,
                               gammastar, gammastar_opt, kmm, kmm_opt, kphio, soilmstress, method_jmaxlim, c4,
                               rd_to_vcmax, beta, c_cost, leafwidth, LAI){
  #---- Fixed parameters--------------------------------------------------------
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  #'
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  
  
  ## viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa)
  ns      <- viscosity_h2o( tcleaf, patm )  # Pa sc4, 1.0,
  ns25    <- viscosity_h2o( kTo, kPo )  # Pa s
  ns_star <- ns / ns25  # (unitless)
  
  # 1. OPTIMAL Vcmax and Jmax----
    
    # Intrinsic quantum efficiency of photosynthesis (phi0)
    phi0_opt = (1/8) *(0.352 + 0.022*tcleaf_opt - 0.00034*tcleaf_opt^(2)) # Temperature dependence function of phi0 (Bernacchi et al.,2003)
    
    # acclimated xiPa (parameter that determines the sensitivity of ci/ca to VPD)
    xiPa = sqrt((beta*(kmm_opt + gammastar_opt))/(1.6*ns_star_opt)) # [Pa^1/2]
    
    # acclimated ci (with acclimated xiPa, and adjusted with the actual VPD)
    ci_opt = (xiPa * ca_opt + gammastar_opt*sqrt(vpd))/
      (xiPa + sqrt(vpd))
    
    
    # OPTIMAL Vcmax
    vcmax_opt  = phi0_opt * ppfd_opt*fapar_opt *
      ((ci_opt + kmm_opt) / (ci_opt + 2*gammastar_opt)) *
      sqrt(1 - (c_cost*(ci_opt + 2*gammastar_opt)/(ci_opt - gammastar_opt))^(2/3)) #[umol m-2 s-1]
    
    # OPTIMAL Jmax
    jmax_opt  = (4 * phi0_opt * ppfd_opt*fapar_opt) / 
      sqrt(1/(1 - (c_cost*( ci_opt + 2*gammastar_opt)/
                     (ci_opt - gammastar_opt))^(2.0/3.0)) - 1) #[umol m-2 s-1]
  
  
  tkleaf_opt = tcleaf_opt + 273.15 # [K]   
  tkleaf = tcleaf + 273.15 # [K]
  
  # 2. INSTANTANEOUS Vcmax and Jmax----
  # The Arrhenius equation constants:
  Ha = 65330# [J mol-1]
  Haj = 43900
  Rgas = 8.314 # [J mol-1 K-1]
  vcmaxAdjusted = vcmax_opt * exp((Ha/Rgas)*(1/tkleaf_opt - 1/tkleaf))
  jmaxAdjusted = jmax_opt * exp((Haj/Rgas)*(1/tkleaf_opt - 1/tkleaf))
  rm(Rgas,Ha,Haj)
  
  # vcmaxAdjusted = vcmax_opt * ftemp_inst_vcmax(tcleaf, tcleaf_opt)
  # jmaxAdjusted = jmax_opt * ftemp_inst_jmax(tcleaf, tcleaf_opt)
  
  
  # 3. instantaneous (with acclimated xiPa, and adjusted with the actual VPD)
  ## photorespiratory compensation point - Gamma-star (Pa)
  gammastar <- calc_gammastar( tcleaf, patm )
  
  # 3.1 Michaelis-Menten coef. (Pa)
  kmm <- calc_kmm( tcleaf, patm )
  
  # 3.2 Intrinsic quantum efficiency of photosynthesis (phi0)
  phi0 = (1/8) *(0.352 + 0.022*tcleaf - 0.00034*tcleaf^(2)) # Temperature dependence function of phi0 (Bernacchi et al.,2003)
  
  # 3.3 ci instantaneous
  ci_inst = (xiPa * ca + gammastar*sqrt(vpd))/(xiPa + sqrt(vpd))
  
  # 3.4 instantaneous chi
  chi_inst = ci_inst/ca
  # ci_inst = chi_inst*ca
  
  # 3.5 Opt chi output
  if(c4){
    out_optchi <- list(
      xi = xiPa,
      chi = chi_inst,
      mc = 1.0,
      mj = 1.0,
      mjoc = 1.0
    )

  } else {
    ## alternative variables
    gamma <- gammastar / ca
    kappa <- kmm / ca

    ## use chi for calculating mj
    mj <- (chi_inst - gamma)/(chi_inst + kappa)

    ## mc
    mc <- (chi_inst - gamma) / (chi_inst +  kappa)

    ## mj:mv
    mjoc <- (chi_inst +  kappa) / (chi_inst + 2.0 * gamma)

    # format output list
    out_optchi <- list(
      xi = xiPa,
      chi = chi_inst,
      mc = mc,
      mj = mj,
      mjoc = mjoc
    )

  }


  # ---- Corrolary preditions ---------------------------------------------------
  # 4.0 intrinsic water use efficiency (in Pa)
  iwue = (ca - ci_inst)/1.6

  #---- Vcmax and light use efficiency -----------------------------------------
  # 7.0 Jmax limitation comes in only at this step
  # if (c4){
  #   out_lue_vcmax <- lue_vcmax_c4(
  #     kphio,
  #     c_molmass,
  #     soilmstress
  #   )
  # 
  # } else if (method_jmaxlim=="wang17"){
  # 
  #   ## apply correction by Jmax limitation
  #   out_lue_vcmax <- lue_vcmax_wang17(
  #     out_optchi,
  #     kphio,
  #     c_molmass,
  #     soilmstress
  #   )
  # 
  # } else if (method_jmaxlim=="smith19"){
  # 
  #   out_lue_vcmax <- lue_vcmax_smith19(
  #     out_optchi,
  #     kphio,
  #     c_molmass,
  #     soilmstress
  #   )
  # 
  # } else if (method_jmaxlim=="none"){
  # 
  #   out_lue_vcmax <- lue_vcmax_none(
  #     out_optchi,
  #     kphio,
  #     c_molmass,
  #     soilmstress
  #   )
  # 
  # } else {
  # 
  #   stop("rpmodel(): argument method_jmaxlim not idetified.")
  # 
  # }
  # 
  
  # 5.0 CALCULATE the assimilation rate

  # acclimated Ac with the acclimated xiPa term
  # a_c = vcmaxAdjusted * out_optchi$mc  #[umol m-2 s-1]
  a_c = vcmaxAdjusted * (ci_inst - gammastar) / (ci_inst + kmm)#[umol m-2 s-1]
  
  # acclimated AJ with the acclimated xiPa term
  # a_j = phi0 * fapar * ppfd * out_optchi$mj * jmaxAdjusted  #[umol m-2 s-1]
  J = (4 *phi0*fapar * ppfd)/sqrt(1 + ((4*phi0*fapar * ppfd)/(jmaxAdjusted))^(2))
  a_j = J/4*(ci_inst - gammastar)/(ci_inst + 2.0 * gammastar)
  
  #---- Corrolary preditions ---------------------------------------------------
  # 6.0 Vcmax25 (vcmax normalized to 25 deg C)
  ftemp25_inst_vcmax  <- ftemp_inst_vcmax( tcleaf, tcleaf_opt, tcref = 25.0 )
  vcmax25  <- vcmax_opt / ftemp25_inst_vcmax
  
  
  ## 7.0  Dark respiration at growth temperature
  ftemp_inst_rd <- ftemp_inst_rd(tcleaf_opt)
  rd  <- rd_to_vcmax * (ftemp_inst_rd / ftemp25_inst_vcmax) * vcmax_opt

  
  # 8.0 Assimilation
  assim <- ifelse(a_j < a_c , a_j, a_c)
  # assim_eq_check <- all.equal(assim, gpp/c_molmass, tol = 0.001)
  # if (! isTRUE(assim_eq_check)) {
  #   warning("rpmodel_subdaily(): Assimilation and GPP are not identical.\n",
  #           assim_eq_check)
  # }
  # Gross Primary Productivity
  gpp <- assim * c_molmass  # in ug C m-2 s-1
  
  ## 9.0  average stomatal conductance
  gs <- assim/(ca - ci_inst)
  e <- 1.6*gs*vpd
  
  ## 10.0 calculate frost and heat cost ####
  #Based in Villar y Merino - 2001 - Comparison of leaf construction costs in woody spe
  # assim <- mapply(calculate_assimilation_lethal, leafwidth, tcleaf, LAI, assim)
  
  ## 11.0 construct list for output
  out <- list(
    gpp             = gpp,   # remove this again later
    assim           = assim,
    ca              = ca,
    gammastar       = gammastar,
    kmm             = kmm,
    chi             = out_optchi$chi,
    xi              = out_optchi$xi,
    mj              = out_optchi$mj,
    mc              = out_optchi$mc,
    ci              = ci_inst,
    iwue            = iwue,
    gs              = gs,
    e               = e,
    vcmax           = vcmaxAdjusted,
    jmax            = jmaxAdjusted,
    vcmax25         = vcmax25,
    rd              = rd,
    tcleaf          = tcleaf,
    vpd_leaf        = vpd,
    ns_star         = ns_star
  )
  
  return(out)
  
  
}


energy_balance <- function(tcleaf_root, tcleaf_opt, vpd_new, ppfd, 
                           ppfd_opt, fapar, fapar_opt, ca, 
                           ca_opt, xi, xiPa, patm, ns_star_opt, 
                           gammastar, gammastar_opt, kmm, 
                           kmm_opt, kphio, soilmstress, 
                           method_jmaxlim, c4, rd_to_vcmax, 
                           beta, c_cost, u, canopy_height,
                           tc, tk, tkleaf, z, LAI, ustar, netrad, mol_gas_const,
                           J_to_mol, lat_heat, mol_mas_wv, sigma, cpm, CP, rho,
                           epsleaf, epssky, frac_PAR, fanir, ei, ea, gb_method, 
                           leafwidth){
  df_res <- rpmodel_jmax_vcmax(tcleaf=tcleaf_root, tcleaf_opt = tcleaf_opt, vpd = vpd_new, ppfd = ppfd, 
                               ppfd_opt = ppfd_opt, fapar = fapar, fapar_opt = fapar_opt, ca = ca, 
                               ca_opt = ca_opt, xi = xi, xiPa = xiPa, patm = patm, ns_star_opt = ns_star_opt, 
                               gammastar = gammastar, gammastar_opt = gammastar_opt, kmm = kmm, 
                               kmm_opt = kmm_opt, kphio = kphio, soilmstress = soilmstress, 
                               method_jmaxlim = method_jmaxlim, c4 = c4, rd_to_vcmax = rd_to_vcmax, 
                               beta = beta, c_cost = c_cost, leafwidth = leafwidth, LAI = LAI)
  
  #Latent Heat Loss calculation
    if(is.na(df_res$gs)){df_res$gs = 0}
    if(length(df_res$gs) == 0){df_res$gs = 0}
    if(is.infinite(df_res$gs)&df_res$gs>0){df_res$gs = 100} 
    if(is.infinite(df_res$gs)&df_res$gs<0){df_res$gs = 0}
    gs = df_res$gs*1.6*1e-6 #stomatal conductance for water
  
    if(is.na(ustar)){
      ust <- calc_ustar(u, canopy_height, z, LAI)
    }else{
      ust <- ustar
    }
  
  # if(!is.na(u)&!is.na(canopy_height)&!is.na(tc)&!is.na(z)&!is.na(LAI)&!is.na(d)){

    gb = calc_ga(ws=u, ust = ust, canopy_height = canopy_height, tcleaf_root,
                 tc, z, LAI, patm,mol_gas_const, rho, tk, gb_method, leafwidth) 
    gb = gb * patm/mol_gas_const/tk #mol m-2 s-1
    gbh = 0.92*gb #boundary layer conductance for heat (Campbell and Norman 1998)
    gbs = gs * gb/(gs + gb)
  
    E = gbs*(vpd_new)*patm/(patm-(ei+ea)/2) #Farquhar and Sharkey 1984e-
    lE = lat_heat*mol_mas_wv*E
  

  #Shortwave Energy Input
    Rs_PAR_Wm2 = fapar*ppfd/(J_to_mol*frac_PAR)
    Rs_NIR_Wm2 = fanir*ppfd/(J_to_mol*frac_PAR) #approximation as for Escobedo et al. 2009 assuming PAR and NIR are equal
    Qsw = Rs_PAR_Wm2 + Rs_NIR_Wm2
    
  #Thermal Infrared Input
    epssky = 1.72 * ((ea*1e-3)/tk)^0.143
    Qtir = epsleaf*epssky*sigma*(tk^4) #sky and air
    
  #Thermal Infra-Red Losses
    Qtirleaf = epsleaf*sigma*tkleaf^4
    # Qtirleaf = 2*epsleaf*sigma*tkleaf^4
    # Qtirleaf = epsleaf*epssky*sigma*tkleaf^4
    
    Rnet = Qsw + Qtir - Qtirleaf

  #Convective Heat Exchange
    Qc = gbh*cpm*(tcleaf_root-tc)
  
  
  return(tibble(Rnet, lE, Qc, Qtir, Qtirleaf, gb, ust))
}



# Define a function that takes leaf width, LAI, assimilation and leaf temperature 
# as inputs and calculates the C cost to replace all the canopy if some temperature 
#thresholds are over-passed
calculate_assimilation_lethal <- function(leafwidth, tcleaf, LAI, assim){

  # Limit lethal temperatures as in Wright et al. 2017; Vitasse et al. 2013 
  # The interaction between freezing tolerance and phenology in temperate deciduous trees;
  #Neuner et al. 2020 Low temperatures at higher elevations require plants to exhibit increased
  #freezing resistance throughout the summer months
  lethal_temp <- (tcleaf < -10 | tcleaf > 50) #LT50

  # Calculate assimilation
  if(lethal_temp){
    
    # Calculate leaf size approximation as from Wright et al. 2017
    leaf_size <- 1.5 * (leafwidth*100)^2 #in cm2
    
    # Calculate g glucose cost per gram of leaf (Villar and Merino 2001 New Phyto)
    # x=c(-2.5,4.483372921615201)
    # y=c(1.7352342158859466,1.2953156822810588)
    # lm(y~x)
    y_g_g <- 1.578 - 0.063 * log(leaf_size)
    
    # Calculate SLA from leaf size
    # x=c(-2.5,4.4654549975989575)
    # y=c(1.6124031007751967,17.240310077519382)
    # lm(y~x)
    y_m2_kg <- 7.221 - 2.244 * log(leaf_size) #in m2 Kg-1
    y_m2_g <- y_m2_kg / 1000 #in m2 g-1
    
    # Convert cost from gram glucose to moles glucose 
    g_glucose_per_to_mole <- 180 #1 mole of glucose weighs 180 g
    mole_glucose_per_mole_C <- 6 #1 mole of glucose contains 6 moles of C
    mole_C_m2 <- y_g_g / y_m2_g / g_glucose_per_to_mole * mole_glucose_per_mole_C
    
    # Calculate total cost
    total_C_cost <- mole_C_m2 * LAI * 1e6 * 0.5 # micromole C Lethal50
    
    assim <- assim - total_C_cost
  }
  
  # Return the assimilation value
  return(assim)
}

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
    TIMESTAMP, tc, vpd, co2, fapar = NA, LAI = NA, ppfd, u = NA, #wind speed in m s^-1
    canopy_height=NA, sw_in = NA, patm = NA, elv = NA, z = NA,
    kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.087182, 0.081785), 0.049977),
    beta = 146.0, c_cost = 0.41, soilm = 1.0, meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.73300,
    c4 = FALSE, method_jmaxlim = "wang17",
    do_ftemp_kphio = TRUE, do_soilmstress = FALSE,do_leaftemp = FALSE,
    energy_params = list(
      epsleaf = 0.96, #thermal absorptivity of the leaf
      ste_bolz = 5.67e-8, #W m^-2 K^-4
      cpm = 75.38, #J mol^-1 ºC-1
      kalb_vis = 0.3, # visible albedo
      kfFEC = 2.0, #Photon flux to energy μmol J-1 (Meek et al., 1984)
      fanir = 0.35 #Fraction of NIR absorbed
    ), returnvar = NULL, verbose = FALSE,
    upscaling_method = c("noon","daily","max_rad"), hour_reference_T = 12, gap_method = "linear",
    xi_acclimated = "on"
){
  
  #---- Fixed parameters--------------------------------------------------------
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  #'
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  
  #parameters
  epsleaf = energy_params["epsleaf"] %>% as.numeric
  sigma = energy_params["ste_bolz"]%>% as.numeric
  cpm = energy_params["cpm"]%>% as.numeric #molar heat capacity of water 
  kalb_vis = energy_params["kalb_vis"]%>% as.numeric
  kfFEC = energy_params["kfFEC"]%>% as.numeric
  fanir = energy_params["fanir"]%>% as.numeric
  lat_heat = 2230 #J g^-1
  mol_gas_const = 8.3144621 #J mol^-1 K^-1
  mol_mas_wv = 18.01528 #g mol-1
  spe_gas_const = mol_gas_const/mol_mas_wv #J g^-1 K^-1
  
  #---- soil moisture stress as a function of soil moisture and mean alpha -----
  if (do_soilmstress) {
    if (length(meanalpha) > 1){
      warning("Argument 'meanalpha' has length > 1. Only the first element is used.")
      meanalpha <- meanalpha[1]
    }
    soilmstress <- calc_soilmstress( soilm, meanalpha, apar_soilm, bpar_soilm )
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
  tibble(TIMESTAMP,tc, vpd, co2, fapar, LAI, ppfd, u, canopy_height, sw_in, patm, meanalpha) %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel(x$tc, x$vpd, x$co2, x$fapar, x$LAI,
                    x$ppfd, x$u, x$canopy_height,x$sw_in, x$patm, 
                    unique(elv), unique(z), kphio, beta,
                    soilm, x$meanalpha, apar_soilm, bpar_soilm, c4,
                    method_jmaxlim, do_ftemp_kphio, do_soilmstress, do_leaftemp = FALSE,
                    energy_params = energy_params)%>% 
      as_tibble()
    return(res)
    }) -> df_Or
  


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
                 z = z,
                 LAI = LAI,
                 canopy_height = canopy_height,
                 patm = patm)
  
  # 2.1 apply the dailyUpscaling function 
  dataDaily <- dailyUpscaling(df = dfIn, 
                                nrWindow = 1, 
                                hour_reference_T = 12, 
                                upscaling_method = "noon")
  
  # 3.0 apply running mean
  dataDaily = runningMean(data_x_running_mean_t = dataDaily, daily_window = 15)
  
  
  # 4.0 Calculate P-model on daily upscaling values
  tibble(dataDaily) %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      res <- rpmodel(x$tc, x$vpd, x$co2, x$fapar, x$LAI,
                     x$ppfd, x$u, x$canopy_height,x$sw_in, x$patm, 
                     unique(elv), unique(z), kphio, beta,
                     soilm, x$meanalpha, apar_soilm, bpar_soilm, c4,
                     method_jmaxlim, do_ftemp_kphio, do_soilmstress, do_leaftemp = FALSE,
                     energy_params = energy_params)%>% 
        as_tibble()
      return(res)
    }) -> df_Mm
  
  
  # 5.0 Downscale from daily to subdaily
  dataDaily = dataDaily %>% cbind(df_Mm) 
  
  DF <- dfIn  %>% 
    cbind(df_Or) %>% 
    left_join(dataDaily %>%
                dplyr::select(-c("YEAR","MONTH","DAY","HOUR","MINUTE")) %>% 
                rename_with( ~ paste0(.x, "_opt"), where(is.numeric)),
              by = "TIMESTAMP")
  
  
  # 5.1 GAP- FILLING
  DF <- lapply(DF, function(y) {
    if (is.numeric(y)) { 
      approx(as.numeric(DF$TIMESTAMP), y, method = gap_method, 
             xout = as.numeric(DF$TIMESTAMP), rule = 2, ties = "constant")$y
      }else{
        y
      }
    }
    ) %>% bind_cols()
  

  
  
  # 6.1 Optimal leaf temperature
  
  tcleaf_new <- DF %>% 
    split(seq(nrow(.))) %>%
    purrr::map_df(function(x){
      es = exp(34.494-4924.99/(x$tc+237.1))/((x$tc+105)^1.57)
      ea = es - x$vpd
      patm = x$patm
      tk = x$tc+273.15
      # tcleaf = x$tc
      # tcleaf_new = tcleaf+1
      tcleaf_new <- uniroot(function(tcleaf_root){
        tkleaf = tcleaf_root+273.15
        # ei = es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tkleaf)) #assuming saturation within the leaf (We don't apply psi_leaf correction es(T)exp(ψleaf V/(RT)))
        ei = exp(34.494-4924.99/(tcleaf_root+237.1))/((tcleaf_root+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
        vpd_new = (ei - ea)
        vpd_new = ifelse(vpd_new<0,0,vpd_new)
        df_res <- rpmodel_jmax_vcmax(tcleaf=tcleaf_root, tcleaf_opt = x$tcleaf_opt, vpd = vpd_new, ppfd = x$ppfd, 
                                     ppfd_opt = x$ppfd_opt, fapar = x$fapar, fapar_opt = x$fapar_opt, ca = x$ca, 
                                     ca_opt = x$ca_opt,  xiPa = x$xiPa, chi_inst = x$chi_inst, ns_star_opt = x$ns_star_opt, 
                                     gammastar = x$gammastar, gammastar_opt = x$gammastar_opt, kmm = x$kmm, 
                                     kmm_opt = x$kmm_opt, kphio = kphio, soilmstress = soilmstress, 
                                     method_jmaxlim = method_jmaxlim, c4 = c4, rd_to_vcmax = rd_to_vcmax, 
                                     xi_acclimated =xi_acclimated,beta = beta, c_cost = c_cost)
        
        #Latent Heat Loss calculation
        if(is.na(df_res$gs)){df_res$gs = 0}
        if(length(df_res$gs) == 0){df_res$gs = 0}
        if(is.infinite(df_res$gs)){df_res$gs = 100} 
        gs = df_res$gs*1.6*1e-6 #stomatal conductance for water
        Hs = gs*cpm*(tcleaf_root-x$tc)
        if(!is.na(x$u)&!is.na(x$canopy_height)&!is.na(x$tc)&!is.na(x$z)&!is.na(x$LAI)){
          # gb = 1/resistance_neutral(ws_mean=x$u, canopy_height = x$canopy_height) * x$patm/mol_gas_const/tk #mol m-2 s-1
          gb = calc_ga(ws=x$u, canopy_height = x$canopy_height,Hs,x$tc,x$z,x$LAI) * x$patm/mol_gas_const/tk #mol m-2 s-1
          gbh = 0.92*gb #boundary layer conductance for heat (Campbell and Norman 1998)
          gbs = gs * gb/(gs + gb)
        }else{
          gbs = gs
          gbh = 0.92*gs
        }
        E = gbs*(vpd_new)*x$patm/(x$patm-(ei+ea)/2) #Farquhar and Sharkey 1984e-
        lE = lat_heat*mol_mas_wv*E
        
        #Shortwave Energy Input
        Rs_PAR_Wm2 = x$ppfd/(kfFEC*(1-kalb_vis))
        Rs_NIR_Wm2 = Rs_PAR_Wm2 #approximation as for Escobedo et al. 2009 assuming PAR and NIR are equal
        Qsw = x$fapar*Rs_PAR_Wm2 + fanir*Rs_NIR_Wm2
        
        #Thermal Infrared Input
        epssky = 1.72 * ((ea*1e-3)/tk)^0.143
        Qtir = epsleaf*epssky*sigma*(tk^4 + tk^4) #sky and air
        
        #Thermal Infra-Red Losses
        # Qtirleaf = epsleaf*sigma*tkleaf^4
        Qtirleaf = 2*epsleaf*sigma*tkleaf^4
        
        #Convective Heat Exchange
        Qc = gbh*cpm*(tcleaf_root-x$tc)
        
        Qsw + Qtir - Qtirleaf - Qc - lE
      },
      c(x$tc-20, x$tc+20))$root
    return(tcleaf_new)
  }) %>% bind_rows() %>% as.numeric()
  tcleaf = tcleaf_new
  tkleaf = tcleaf+273.15
  es = exp(34.494-4924.99/(tc+237.1))/((tc+105)^1.57)
  ea = es - vpd
  ei = exp(34.494-4924.99/(tcleaf_new+237.1))/((tcleaf_new+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
  vpd_new = (ei - ea)
  vpd_new = ifelse(vpd_new<0,0,vpd_new)

  out <- rpmodel_jmax_vcmax(tcleaf=tcleaf_new, tcleaf_opt = DF$tcleaf_opt, vpd = vpd_new, ppfd = DF$ppfd, ppfd_opt = DF$ppfd_opt,
                            fapar = DF$fapar, fapar_opt = DF$fapar_opt, ca = DF$ca, ca_opt = DF$ca_opt, 
                            xiPa = DF$xiPa, chi_inst = DF$chi_inst, ns_star_opt = DF$ns_star_opt, 
                            gammastar = DF$gammastar, gammastar_opt = DF$gammastar_opt, kmm = DF$kmm, kmm_opt = DF$kmm_opt,
                            kphio = kphio, soilmstress = soilmstress, method_jmaxlim = method_jmaxlim, c4 = c4, rd_to_vcmax = rd_to_vcmax,
                            xi_acclimated =xi_acclimated, beta = beta, c_cost = c_cost)
  
  return(out)

}






dailyUpscaling <- function(df = dfIn, nrWindow = 1, hour_reference_T = 12, 
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
  if (upscaling_method == "max_rad")
    colMandatory = c(colMandatory,"sw_in")
  
  if (!headerControl_dd(df = df, colMandatory = colMandatory))
    stop(headerControl_dd(df = df, colMandatory = colMandatory, showMsg = TRUE))
  
  #2.0 T calculation
  df <- df %>% as.data.frame()
  dfDayT = df[1,]    
  
  
  for (cicloAnno in sort(unique(df$YEAR))) {
    for (cicloMesi in seq(1,12)) {
      for (cicloGiorni in seq(1,31)) {
        
        posDay = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi &
                         df$DAY == cicloGiorni)
        if (length(posDay) == 0) next
        
        for (hourReference in hour_reference_T) {
          # upscaling_method noon
          if (upscaling_method == "noon") {
            posReference = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                   df$HOUR == hourReference & df$MINUTE == 0)
            if (length(posReference) == 0) next
            windowDay = seq(posReference - nrWindow,posReference + nrWindow)
          }
          # upscaling_method daily
          if (upscaling_method == "daily") {
            posReference = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                   df$HOUR == 12 & df$MINUTE == 0)
            windowDay = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni)
          }
          # upscaling_method max_rad
          if (upscaling_method == "max_rad") {
            posReferenceMax = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
                                      is.na(df$SWINPOT) == 0)
            if ( length(posReferenceMax) > 0 ) {
              
              posMaxSWINPOT = which(df$SWINPOT[posReferenceMax] == max(df$SWINPOT[posReferenceMax],na.rm = T)[1])
              
              posReference = posReferenceMax[posMaxSWINPOT[1]]
              windowDay = seq(posReference - nrWindow,posReference + nrWindow)
            } else {
              posReference = which(df$YEAR == cicloAnno & df$MONTH == cicloMesi & df$DAY == cicloGiorni &
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
                dfDay[1,col] = mean(tmp)
              }
              rm(tmp,posNa)
            }else{next}
          }
          
          dfDay$TIMESTAMP = df$TIMESTAMP[posReference]
          
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


runningMean <- function(data_x_running_mean_t = df, daily_window = 10) {
  
  unique_hour = sort(unique(data_x_running_mean_t$HOUR))
  
  data_running_mean_t = data_x_running_mean_t[1,]
  
  for (cy_hour_ref in unique_hour) {
    pos_hour = which(data_x_running_mean_t$HOUR == cy_hour_ref)
    
    data_x_running_mean = data_x_running_mean_t[pos_hour,]
    rm(pos_hour)
    # window of time 
    lunghezza_finestra = daily_window - 1# remotion of one day since the computation starts from the 1st
    data_running_mean = data_x_running_mean[1,]
    
    for ( ciclo_inizio_mm in seq(2,nrow(data_x_running_mean)) ) {
      # definition of positions to compute the running mean by making a counter in reverse
      posizioni_rm = (ciclo_inizio_mm - lunghezza_finestra):ciclo_inizio_mm
      
      pos_meno1 = which(posizioni_rm < 1)
      if ( length(pos_meno1) > 0 ) posizioni_rm = posizioni_rm[-1*pos_meno1]
      rm(pos_meno1)
      
      # cycle for computing the mean of the dataset's variables
      media_1 = data_x_running_mean[1,]
      for ( ciclo_variabili in colnames(media_1) ) {
        # dataset creation with values to use
        data_1 = data_x_running_mean[posizioni_rm,ciclo_variabili]
        if ( !is.numeric(data_1) ) {
          media_1[ciclo_variabili] = NA
          next
        }
        # NA remotion
        pos_na = which(is.na(data_1) == 1 )
        if ( length(pos_na) > 0 ) data_1 = data_1[-1*pos_na]
        rm(pos_na)
        # if there are missing values it puts NA, otherwise it computes the mean
        if ( length(data_1) == 0 ) {
          media_1[ciclo_variabili] = NA
        } else {
          media_1[ciclo_variabili] = mean(data_1)
        }
        rm(data_1)
      }
      rm(ciclo_variabili)
      media_1$TIMESTAMP = data_x_running_mean$TIMESTAMP[ciclo_inizio_mm]
      data_running_mean = rbind(data_running_mean,media_1)
      rm(media_1)
      rm(posizioni_rm)
    }
    rm(ciclo_inizio_mm)
    rm(lunghezza_finestra)
    rm(data_x_running_mean)
    data_running_mean_t = rbind(data_running_mean_t,data_running_mean)
  }
  
  data_running_mean_t = data_running_mean_t[-1,]
  
  # library(lubridate)
  # data_running_mean_t$YEAR = year(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
  # data_running_mean_t$MONTH = month(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
  # data_running_mean_t$DAY = day(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
  # data_running_mean_t$HOUR = hour(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
  # data_running_mean_t$MINUTE = minute(ymd_hm(as.character(data_running_mean_t$TIMESTAMP)))
  
  data_running_mean_t = data_running_mean_t[order(data_running_mean_t$TIMESTAMP),]
  
  return(data_running_mean_t)
}



rpmodel_jmax_vcmax <- function(tcleaf, tcleaf_opt, vpd, ppfd, ppfd_opt, fapar, fapar_opt, ca, ca_opt, xiPa, chi_inst, ns_star_opt,
                               gammastar, gammastar_opt, kmm, kmm_opt, kphio, soilmstress, method_jmaxlim, c4,
                               rd_to_vcmax, xi_acclimated, beta, c_cost){
  #---- Fixed parameters--------------------------------------------------------
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  #'
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  
  
  # 5.2 OPTIMAL Vcmax and Jmax----
  # 5.2.1 Vcmax with acclimated 'xiPa', 'ci' and 'phi0' (ON)
  if (xi_acclimated == 'off') {
    cat(sprintf('vcmax_opt and jmax_opt from upscaling\n'))
    DF[,'ci'] = NA
    DF[,'xiPa'] = NA
  } else {
    
    # cat(sprintf('vcmax_opt and jmax_opt calculated\n'))
    
    # Intrinsic quantum efficiency of photosynthesis (phi0)
    phi0 = (1/8) *(0.352 + 0.022*tcleaf_opt - 0.00034*tcleaf_opt^(2)) # Temperature dependence function of phi0 (Bernacchi et al.,2003)
    
    # acclimated xiPa (parameter that determines the sensitivity of ci/ca to VPD)
    xiPa = sqrt((beta*(kmm_opt + gammastar_opt))/(1.6*ns_star_opt)) # [Pa^1/2]
    
    # acclimated ci (with acclimated xiPa, and adjusted with the actual VPD)
    ci = (xiPa * ca_opt + gammastar_opt*sqrt(vpd))/
      (xiPa + sqrt(vpd))
    
    
    # OPTIMAL Vcmax
    vcmax_opt  = phi0 * ppfd_opt*fapar_opt *
      ((ci + kmm_opt) / (ci + 2*gammastar_opt)) *
      sqrt(1 - (c_cost*(ci + 2*gammastar_opt)/(ci - gammastar_opt))^(2/3)) #[umol m-2 s-1]
    
    # OPTIMAL Jmax
    jmax_opt  = (4 * phi0 * ppfd_opt*fapar_opt) / 
      sqrt(1/(1 - (c_cost*( ci + 2*gammastar_opt)/
                     (ci - gammastar_opt))^(2.0/3.0)) - 1) #[umol m-2 s-1]
  }
  
  
  tkleaf_opt = tcleaf_opt + 273.15 # [K]   
  tkleaf = tcleaf + 273.15 # [K]
  
  # 5.2.3 INSTANTANEOUS Vcmax and Jmax----
  # The Arrhenius equation constants:
  Ha = 65330# [J mol-1]
  Haj = 43900
  Rgas = 8.314 # [J mol-1 K-1]
  
  
  vcmaxAdjusted = vcmax_opt * exp((Ha/Rgas)*(1/tkleaf_opt - 1/tkleaf))
  jmaxAdjusted = jmax_opt * exp((Haj/Rgas)*(1/tkleaf_opt - 1/tkleaf))
  
  rm(Rgas,Ha,Haj)
  
  
  # 5.3 instantaneous ci (with acclimated xiPa, and adjusted with the actual VPD)
  ci_inst = (xiPa * ca + gammastar*sqrt(vpd))/(xiPa + sqrt(vpd))
  
  # 5.4 instantaneous chi
  chi_inst = ci_inst/ca 
  
  # 5.5 Opt chi output
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
  # 6.0 intrinsic water use efficiency (in Pa)
  iwue = (ca - ci_inst)/1.6
  
  #---- Vcmax and light use efficiency -----------------------------------------
  # 7.0 Jmax limitation comes in only at this step
  if (c4){
    out_lue_vcmax <- lue_vcmax_c4(
      kphio,
      c_molmass,
      soilmstress
    )
    
  } else if (method_jmaxlim=="wang17"){
    
    ## apply correction by Jmax limitation
    out_lue_vcmax <- lue_vcmax_wang17(
      out_optchi,
      kphio,
      c_molmass,
      soilmstress
    )
    
  } else if (method_jmaxlim=="smith19"){
    
    out_lue_vcmax <- lue_vcmax_smith19(
      out_optchi,
      kphio,
      c_molmass,
      soilmstress
    )
    
  } else if (method_jmaxlim=="none"){
    
    out_lue_vcmax <- lue_vcmax_none(
      out_optchi,
      kphio,
      c_molmass,
      soilmstress
    )
    
  } else {
    
    stop("rpmodel(): argument method_jmaxlim not idetified.")
    
  }
  
  
  # 8.0 CALCULATE the assimilation rate
  # acclimated Ac with the acclimated xiPa term
  a_c = vcmaxAdjusted * out_optchi$mc  #[umol m-2 s-1]
  
  # acclimated AJ with the acclimated xiPa term
  a_j = kphio * fapar*ppfd * out_optchi$mj * jmaxAdjusted  #[umol m-2 s-1]
  
  
  
  #---- Corrolary preditions ---------------------------------------------------
  # 9.0 Vcmax25 (vcmax normalized to 25 deg C)
  ftemp25_inst_vcmax  <- ftemp_inst_vcmax( tcleaf, tcleaf_opt, tcref = 25.0 )
  vcmax25_unitiabs  <- out_lue_vcmax$vcmax_unitiabs / ftemp25_inst_vcmax
  
  
  ## 10.0  Dark respiration at growth temperature
  ftemp_inst_rd <- ftemp_inst_rd(tcleaf_opt)
  rd_unitiabs  <- rd_to_vcmax * (ftemp_inst_rd / ftemp25_inst_vcmax) * out_lue_vcmax$vcmax_unitiabs
  
  ## Dark respiration
  rd <- ppfd * fapar * rd_unitiabs
  
  # 11.0 Assimilation
  assim <- ifelse(a_j < a_c , a_j, a_c)
  # assim_eq_check <- all.equal(assim, gpp/c_molmass, tol = 0.001)
  # if (! isTRUE(assim_eq_check)) {
  #   warning("rpmodel_subdaily(): Assimilation and GPP are not identical.\n",
  #           assim_eq_check)
  # }
  # Gross Primary Productivity
  gpp <- assim * c_molmass  # in ug C m-2 s-1
  
  ## 12.0  average stomatal conductance
  gs <- assim/(ca - ci_inst)
  e <- 1.6*gs*vpd
  
  ## 13.0 construct list for output
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
    rd              = rd,
    tcleaf          = tcleaf,
    vpd_leaf        = vpd
  )
  
  return(out)
  
  
}


# Function to fill in the missing data of a vector according to the settings
# param: v, data vector to be filled
# param: vYear, vector of the years
# param: vMonth, vector of the months
# param: vDay, vector of the days
# param: vHour, vector of the hours
# param: vMinute, vector of the minutes
# param: approccio (approach) if 'constant' replicates measured data on missing values;
#                             if 'linear' it applies linear regression to fill in missing values among measured data
# param: showMsg (if T shows the function messages)
# return: a vector with required inputs for pmodelPlus
# rdname: gapFilling

# gapFilling = function(v = vettore,v_name = NA,
#                       vYear = vYear, vMonth = vMonth, vDay = vDay,
#                       vHour = vHour, vMinute = vMinute, showMsg = F,
#                       approccio = 'constant') {
#   
#   # check if there are values in v 
#   if ( sum(is.na(v)) == length(v) ) {
#     warning('variable not gap filled\n') 
#     return(v)
#   }
#   str_msg = ''
#   if (!is.na(v_name))
#     str_msg = paste0(str_msg,sprintf('GAP FILLING variable: %s\n',v_name))
#   
#   str_msg = paste0(str_msg,
#                    sprintf('approach: %s\nmissing value%s: %s\n',
#                            approccio,
#                            ifelse(sum(is.na(v)) == 1,'','s'),sum(is.na(v))))
#   
#   if (showMsg)
#     cat(str_msg)
#   
#   if (approccio == 'constant') {
#     # 
#     # for (cyr in seq(2,length(v))) {
#     # 
#     #   if (is.na(v[cyr]))
#     #     v[cyr] = v[cyr-1]
#     # }
#     # rm(cyr)
#     posValori = which(is.na(v) == 0)
#     for (contaV in posValori) {
#       
#       yearToUse  = vYear[contaV]
#       monthToUse = vMonth[contaV]
#       dayToUse   = vDay[contaV]
#       vToUse = v[contaV]
#       
#       posToUse = which(
#         vYear == yearToUse &
#           vMonth == monthToUse &
#           vDay == dayToUse)
#       
#       if ( length(posToUse) > 0 )
#         v[posToUse] = vToUse
#       
#       rm(yearToUse,monthToUse,dayToUse)
#       rm(vToUse,posToUse)
#     }
#     rm(contaV)
#   }
#   #approccio: 'linear' 
#   if (approccio == 'linear') {
#     posValori = which(is.na(v) == 0)
#     for (contaV in seq(2,length(posValori))) {
#       # slope computetation 
#       x1 = posValori[contaV - 1]
#       x2 = posValori[contaV]
#       y1 = v[x1]
#       y2 = v[x2]
#       slope = (y2 - y1) / (x2-x1)
#       intercept = y1 - (slope*x1)
#       itp = (slope * seq(x1,x2)) + intercept
#       v[seq(x1,x2)] = itp
#     }
#   }
#   
#   if (showMsg)
#     cat(sprintf('GAP FILLING COMPLETED\nmissing value%s:%d\n',
#                 ifelse(sum(is.na(v)) == 1,'','s'),sum(is.na(v))))
#   
#   return(v)
#   
# }

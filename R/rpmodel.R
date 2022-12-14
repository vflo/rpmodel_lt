



rpmodel <- function(
    tc,
    vpd,
    co2,
    fapar,
    ppfd,
    u = NA,
    canopy_height = NA,
    patm = NA,
    elv = NA,
    kphio = ifelse(c4, 1.0,
                   ifelse(do_ftemp_kphio,
                          ifelse(do_soilmstress,
                                 0.087182,
                                 0.081785),
                          0.049977)),
    beta = ifelse(c4, 146/9, 146),
    soilm = stopifnot(!do_soilmstress),
    meanalpha = 1.0,
    apar_soilm = 0.0,
    bpar_soilm = 0.73300,
    c4 = FALSE,
    method_jmaxlim = "wang17",
    do_ftemp_kphio = TRUE,
    do_soilmstress = FALSE,
    do_leaftemp = FALSE,
    energy_params = list(
      epsleaf = 0.96, #thermal absorptivity of the leaf
      ste_bolz = 5.67e-8, #W m^-2 K^-4
      cpm = 75.38, #J mol^-1 ºC-1
      kalb_vis = 0.3, # visible albedo
      kfFEC = 2.0, #Photon flux to energy μmol J-1 (Meek et al., 1984)
      fanir = 0.35 #Fraction of NIR absorbed
    ),
    returnvar = NULL,
    verbose = FALSE){
  
  if(do_leaftemp == FALSE){
    #set tcleaf equal to ambient growth temperature
    tcleaf = tc
    
    rpmodel_core(
      tc,
      tcleaf,
      vpd,
      co2,
      fapar,
      ppfd,
      patm = patm,
      elv = elv,
      kphio = ifelse(c4, 1.0,
                     ifelse(do_ftemp_kphio,
                            ifelse(do_soilmstress,
                                   0.087182,
                                   0.081785),
                            0.049977)),
      beta = ifelse(c4, 146/9, 146),
      soilm = stopifnot(!do_soilmstress),
      meanalpha = 1.0,
      apar_soilm = 0.0,
      bpar_soilm = 0.73300,
      c4 = FALSE,
      method_jmaxlim = "wang17",
      do_ftemp_kphio = TRUE,
      do_soilmstress = FALSE,
      returnvar = NULL,
      verbose = FALSE)
    
  }else{
    
    rpmodel_lt(
      tc,
      vpd,
      co2,
      fapar,
      ppfd,
      u, #wind speed in m s^-1
      canopy_height,
      patm = patm,
      elv = elv,
      kphio = ifelse(c4, 1.0,
                     ifelse(do_ftemp_kphio,
                            ifelse(do_soilmstress,
                                   0.087182,
                                   0.081785),
                            0.049977)),
      beta = ifelse(c4, 146/9, 146),
      soilm = stopifnot(!do_soilmstress),
      meanalpha = 1.0,
      apar_soilm = 0.0,
      bpar_soilm = 0.73300,
      c4 = FALSE,
      method_jmaxlim = "wang17",
      do_ftemp_kphio = TRUE,
      do_soilmstress = FALSE,
      do_leaftemp = FALSE,
      energy_params = energy_params,
      returnvar = NULL,
      verbose = FALSE)
    
    
  }
  
  
}




rpmodel_lt <- function(
    tc,
    vpd,
    co2,
    fapar,
    ppfd,
    u, #wind speed in m s^-1
    canopy_height,
    patm = NA,
    elv = NA,
    kphio = ifelse(c4, 1.0,
                   ifelse(do_ftemp_kphio,
                          ifelse(do_soilmstress,
                                 0.087182,
                                 0.081785),
                          0.049977)),
    beta = ifelse(c4, 146/9, 146),
    soilm = stopifnot(!do_soilmstress),
    meanalpha = 1.0,
    apar_soilm = 0.0,
    bpar_soilm = 0.73300,
    c4 = FALSE,
    method_jmaxlim = "wang17",
    do_ftemp_kphio = TRUE,
    do_soilmstress = FALSE,
    do_leaftemp = FALSE,
    energy_params = list(
      epsleaf = 0.96, #thermal absorptivity of the leaf
      ste_bolz = 5.67e-8, #W m^-2 K^-4
      cpm = 75.38, #J mol^-1 ºC-1
      kalb_vis = 0.3, # visible albedo
      kfFEC = 2.0, #Photon flux to energy μmol J-1 (Meek et al., 1984)
      fanir = 0.35 #Fraction of NIR absorbed
    ),
    returnvar = NULL,
    verbose = FALSE){
  
  #parameters
  epsleaf = energy_params["epsleaf"] %>% as.numeric
  sigma = energy_params["ste_bolz"]%>% as.numeric
  cpm = energy_params["cpm"]%>% as.numeric #molar heat capacity of water 
  kalb_vis = energy_params["kalb_vis"]%>% as.numeric
  kfFEC = energy_params["kfFEC"]%>% as.numeric
  fanir = energy_params["fanir"]%>% as.numeric
  tk = tc+273.15
  lat_heat = 2230 #J g^-1
  mol_gas_const = 8.3144621 #J mol^-1 K^-1
  mol_mas_wv = 18.01528 #g mol-1
  spe_gas_const = mol_gas_const/mol_mas_wv #J g^-1 K^-1
  es_T0 = 610.8 #Pa
  # es = es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tk)) # Pa (Clausius–Clapeyron relation) 
  es = exp(34.494-4924.99/(tc+237.1))/((tc+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
  ea = es - vpd
  # Check arguments
  if (identical(NA, elv) && identical(NA, patm)){
    stop(
      "Aborted. Provide either elevation (arugment elv) or 
     atmospheric pressure (argument patm)."
    )
  } else if (!identical(NA, elv) && identical(NA, patm)){
    if (verbose) {
      warning(
        "Atmospheric pressure (patm) not provided. Calculating it as a
      function of elevation (elv), assuming standard atmosphere 
      (101325 Pa at sea level)."
      )
    }
    
    patm <- calc_patm(elv)
  }
  
  #initialize recursive 
  tcleaf = tc
  tcleaf_new = tcleaf+1
  
  # while(abs(tcleaf-tcleaf_new)>1e-3){
    # tcleaf = tcleaf_new
    # tkleaf = tcleaf+273.15
    # vpd_new = (es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tkleaf)) - ea)
  #   while(vpd_new < 0) {
  #     tcleaf_new = (tcleaf_new+tc)/2 + runif(1)
  #     tcleaf = tcleaf_new
  #     tkleaf = tcleaf+273.15
  #     vpd_new = (es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tkleaf)) - ea)
  #     # vpd_new = 0.01
  #   }
    # calculate pmodel

    
    tcleaf_new <- uniroot(function(tcleaf_root){
      tkleaf = tcleaf_root+273.15
      # ei = es_T0*exp(lat_heat/spe_gas_const*(1/273.15-1/tkleaf)) #assuming saturation within the leaf (We don't apply psi_leaf correction es(T)exp(ψleaf V/(RT)))
      ei = exp(34.494-4924.99/(tcleaf_root+237.1))/((tcleaf_root+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
      vpd_new = (ei - ea)
      vpd_new = ifelse(vpd_new<0,0,vpd_new)
      df_res <- rpmodel_core(
        tc,
        tcleaf_root,
        vpd_new,
        co2,
        fapar,
        ppfd,
        patm = patm,
        elv = elv,
        kphio = ifelse(c4, 1.0,
                       ifelse(do_ftemp_kphio,
                              ifelse(do_soilmstress,
                                     0.087182,
                                     0.081785),
                              0.049977)),
        beta = ifelse(c4, 146/9, 146),
        soilm = stopifnot(!do_soilmstress),
        meanalpha = 1.0,
        apar_soilm = 0.0,
        bpar_soilm = 0.73300,
        c4 = FALSE,
        method_jmaxlim = "wang17",
        do_ftemp_kphio = TRUE,
        do_soilmstress = FALSE,
        returnvar = NULL,
        verbose = FALSE)
      #Latent Heat Loss calculation
      if(is.na(df_res$gs)){df_res$gs = 0}
      if(length(df_res$gs) == 0){df_res$gs = 0}
      if(is.infinite(df_res$gs)){df_res$gs = 100} 
      gs = df_res$gs*1.6*1e-6 #stomatal conductance for water
      if(!is.na(u)&!is.na(canopy_height)){
        gb = 1/resistance_neutral(ws_mean=u, canopy_height = canopy_height) *patm/mol_gas_const/tk #mol m-2 s-1
        gbh = 0.92*gb #boundary layer conductance for heat (Campbell and Norman 1998)
        gbs = gs * gb/(gs + gb)
      }else{
        gbs = gs
        gbh = 0.92*gs
      }
      E = gbs*(vpd_new)*patm/(patm-(ei+ea)/2) #Farquhar and Sharkey 1984e-
      lE = lat_heat*mol_mas_wv*E
      
      #Shortwave Energy Input
      Rs_PAR_Wm2 = ppfd/(kfFEC*(1-kalb_vis))
      Rs_NIR_Wm2 = Rs_PAR_Wm2 #approximation as for Escobedo et al. 2009 assuming PAR and NIR are equal
      Qsw = fapar*Rs_PAR_Wm2 + fanir*Rs_NIR_Wm2
      
      #Thermal Infrared Input
      epssky = 1.72 * ((ea*1e-3)/tk)^0.143
      Qtir = epsleaf*epssky*sigma*tk^4
      
      #Thermal Infra-Red Losses
      Qtirleaf = epsleaf*sigma*tkleaf^4
      # Qtirleaf = 2*epsleaf*sigma*tkleaf^4
      
      #Convective Heat Exchange
      Qc = gbh*cpm*(tcleaf_root-tc)
      
      Qsw + Qtir - Qtirleaf - Qc - lE
    },
    c(tc-15, tc+15))$root
    
    tcleaf = tcleaf_new
    tkleaf = tcleaf+273.15
    ei = exp(34.494-4924.99/(tcleaf_new+237.1))/((tcleaf_new+105)^1.57) #Pa https://journals.ametsoc.org/view/journals/apme/57/6/jamc-d-17-0334.1.xml
    vpd_new = (ei - ea)
    df_res <- rpmodel_core(
      tc,
      tcleaf_new,
      vpd_new,
      co2,
      fapar,
      ppfd,
      patm = patm,
      elv = elv,
      kphio = ifelse(c4, 1.0,
                     ifelse(do_ftemp_kphio,
                            ifelse(do_soilmstress,
                                   0.087182,
                                   0.081785),
                            0.049977)),
      beta = ifelse(c4, 146/9, 146),
      soilm = stopifnot(!do_soilmstress),
      meanalpha = 1.0,
      apar_soilm = 0.0,
      bpar_soilm = 0.73300,
      c4 = FALSE,
      method_jmaxlim = "wang17",
      do_ftemp_kphio = TRUE,
      do_soilmstress = FALSE,
      returnvar = NULL,
      verbose = FALSE)

  # }
  
  #return the result
  return(df_res)
  
}

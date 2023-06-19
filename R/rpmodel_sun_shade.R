rpmodel_sun_shade <- function(
    tc,
    tcleaf = NA,
    vpd,
    co2,
    fapar,
    ppfd_sun,
    ppfd_shade,
    patm = NA,
    elv = NA,
    beta = ifelse(c4, 146/9, 146),
    c_cost = 0.41, # Change to 0.05336251 if method_jmalim = "smith19" is used
    soilm = stopifnot(!do_soilmstress),
    AI = 1,
    c4 = FALSE,
    method_jmaxlim = "wang17",
    do_ftemp_kphio = TRUE,
    do_soilmstress = FALSE,
    returnvar = NULL,
    verbose = FALSE 
){
  res_sun <- rpmodel_core(tc,
               tcleaf,
               vpd,
               co2,
               1,
               ppfd_sun,
               patm,
               elv,
               kphio=ifelse(c4, 1.0,
                            ifelse(do_ftemp_kphio,
                                   ifelse(do_soilmstress,
                                          0.087182,
                                          0.081785),
                                   0.049977)),
               beta = ifelse(c4, 146/9, 146),
               c_cost, # Change to 0.05336251 if method_jmalim = "smith19" is used
               soilm,
               AI ,
               c4 ,
               method_jmaxlim,
               do_ftemp_kphio,
               do_soilmstress,
               returnvar,
               verbose )
  
  
  res_shade <- rpmodel_core(tc,
                          tcleaf,
                          vpd,
                          co2,
                          1,
                          ppfd_shade,
                          patm,
                          elv,
                          kphio=ifelse(c4, 1.0,
                                       ifelse(do_ftemp_kphio,
                                              ifelse(do_soilmstress,
                                                     0.087182,
                                                     0.081785),
                                              0.049977)),
                          beta = ifelse(c4, 146/9, 146),
                          c_cost, # Change to 0.05336251 if method_jmalim = "smith19" is used
                          soilm,
                          AI ,
                          c4 ,
                          method_jmaxlim,
                          do_ftemp_kphio,
                          do_soilmstress,
                          returnvar,
                          verbose )
  
  res_sun_names <- paste0(names(res_sun),"_sun")
  names(res_sun) <- res_sun_names
  res_shade_names<- paste0(names(res_shade),"_shade")
  names(res_shade) <- res_shade_names
  res <- c(res_sun,res_shade)
  return(res)
}

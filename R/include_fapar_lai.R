
include_fapar_lai <- function(fn, fapar_noaa, fapar){
  fapar_noaa_filled <- fapar_noaa %>% 
    filter(site == unique(fn$site),
           FparLai_QC <= 1
    )%>%
    group_by(date,site) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(timestamp = date, Fapar = Fpar,Lai = Lai) %>%
    mutate(Fapar = Fapar*0.001, Lai = Lai*0.001, source = "NOAA")
  
  fapar_filled <- fapar %>% 
    filter(site == unique(fn$site),
           FparLai_QC %in% c(0,2)
    ) %>%
    group_by(date,site) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(timestamp = date, Fapar = Fpar,Lai = Lai) %>%
    mutate(Fapar = Fapar*0.01,Lai = Lai*0.1, source = "MODIS")
  
  fapar_filled <- fapar_noaa_filled %>% bind_rows(fapar_filled) %>% 
    filter(timestamp>= min(fn$timestamp, na.rm = TRUE),timestamp< max(fn$timestamp, na.rm = TRUE))
  
  min_date = min(fapar_filled$timestamp %>% 
                   lubridate::as_date())
  
  # Build FAPAR model
  k_ <- as.integer(nrow(fapar_filled)*0.2)
  if(k_ >= 3){
    # x <- fapar_filled %>%
    #   dplyr::select(timestamp,Fapar) %>%
    #   mutate(timestamp = timestamp %>% lubridate::ymd() %>% 
    #            as.numeric()) %>% select(timestamp)
    # model_fapar <- splinefun(x=x$timestamp,y=fapar_filled$Fapar,method ='natural')
    # model_fapar <- qgam::qgam(Fapar ~ s(timestamp, bs = 'tp',k = k_), data = fapar_filled %>%
    #                             dplyr::select(timestamp,Fapar) %>%
    #                             mutate(timestamp = timestamp %>% lubridate::ymd() %>%
    #                                      as.numeric()), qu = 0.8)
    # 
    # model_lai <- qgam::qgam(Lai ~ s(timestamp, bs = 'tp', k = k_), data = fapar_filled %>%
    #                           dplyr::select(timestamp,Lai) %>%
    #                           mutate(timestamp = timestamp %>% lubridate::ymd() %>%
    #                                    as.numeric()), qu = 0.8)
    # x <- fapar_filled$timestamp %>%lubridate::ymd() %>% as.numeric()
    # y_fapar <- fapar_filled$Fapar
    # y_lai <- fapar_filled$Lai
    # model_fapar <- GauPro::GauPro(x, y_fapar, parallel=TRUE)
    # model_lai <- GauPro::GauPro(x, y_lai, parallel=TRUE)
    # # 
    # pred_FAPAR <- model_fapar$pred(fn[,"timestamp"] %>%
    #                                  lubridate::as_date() %>%
    #                                  as.numeric())
    # pred_LAI <- model_lai$pred(fn[,"timestamp"] %>%
    #                              lubridate::as_date() %>%
    #                              as.numeric())
    # pred_FAPAR <- predict(model_fapar, tibble(timestamp = fn$timestamp %>% 
    #                                             lubridate::as_date() %>%
    #                                             as.numeric()%>% 
    #                                             unique()))
    # pred_LAI <- predict(model_lai, tibble(timestamp = fn$timestamp %>%
    #                                         lubridate::as_date() %>%
    #                                         as.numeric()%>% 
    #                                         unique()))
    
    z = forecastML::fill_gaps(fapar_filled %>%
                                dplyr::select(timestamp,Fapar,Lai) %>%
                                mutate(date = timestamp %>% 
                                         lubridate::ymd()) %>% 
                                ungroup() %>% 
                                dplyr::select(date,Fapar,Lai), 
                              date_col = 1,
                              frequency = '1 day')
    z <- z %>% mutate(Fapar = as.numeric(Fapar),
                      Lai = as.numeric(Lai))
    FAPAR = zoo::na.approx(z[,c(2)],na.rm=FALSE, maxgap = as.numeric(100))
    LAI = zoo::na.approx(z[,c(3)],na.rm=FALSE, maxgap = as.numeric(100))
    FAPAR = ifelse(FAPAR<0,0,FAPAR)
    FAPAR = ifelse(FAPAR>1,1,FAPAR)
    LAI = ifelse(LAI<0,0,LAI)
    z <- z %>% bind_cols(FAPAR= FAPAR, LAI = LAI)
    
    fn %>% 
      mutate(date = lubridate::as_date(timestamp)) %>% 
      left_join(z %>% dplyr::select(-c(Fapar, Lai)),
      # left_join(tibble(date = fn$timestamp%>% 
      #                    lubridate::as_date() %>% 
      #                    unique(), 
      #                  FAPAR = FAPAR) %>% 
                  # filter(date > min_date) %>%
                  # mutate(FAPAR = case_when(FAPAR<0~0,
                  #                          FAPAR>1~1,
                  #                          TRUE~FAPAR)),
                by = "date") %>% 
      # left_join(tibble(date = fn$timestamp%>% 
      #                    lubridate::as_date() %>% 
      #                    unique(),  
      #                  LAI = LAI) %>%
      #             filter(date > min_date) %>% 
      #             mutate(LAI = case_when(LAI<0~0,
      #                                    LAI>=0~LAI)),
      #           by = "date") %>% 
      left_join(fapar_filled %>% 
                  dplyr::mutate(date = lubridate::ymd(timestamp)) %>% 
                  dplyr::ungroup() %>% 
                  dplyr::select(-timestamp), 
                by = "date")
    
  }else{
    fn %>%
      left_join(fapar_filled %>% mutate(timestamp = lubridate::ymd(timestamp)), by = "timestamp") %>% 
      mutate(FAPAR = mean(Fapar,na.rm = TRUE),
             LAI = mean(Lai,na.rm = TRUE))
  }
}



include_fapar_lai <- function(sfn, fapar_noaa, fapar){
  fapar_noaa_filled <- fapar_noaa %>% 
    filter(si_code == unique(sfn$si_code),
           FparLai_QC <= 1
    )%>%
    group_by(date,si_code) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(timestamp = date, Fapar = Fpar,Lai = Lai) %>%
    mutate(Fapar = Fapar*0.001, Lai = Lai*0.001, source = "NOAA")
  
  fapar_filled <- fapar %>% 
    filter(si_code == unique(sfn$si_code),
           FparLai_QC %in% c(0,2)
    ) %>%
    group_by(date,si_code) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(timestamp = date, Fapar = Fpar,Lai = Lai) %>%
    mutate(Fapar = Fapar*0.01,Lai = Lai*0.1, source = "MODIS")
  
  fapar_filled <- fapar_noaa_filled %>% bind_rows(fapar_filled) %>% 
    filter(timestamp>= min(sfn$timestamp, na.rm = TRUE),timestamp< max(sfn$timestamp, na.rm = TRUE))
  
  min_data = min(fapar_filled$timestamp)
  
  # Build FAPAR model
  k_ <- as.integer(nrow(fapar_filled)*0.2)
  if(k_ >= 3){
    # x <- fapar_filled %>%
    #   dplyr::select(timestamp,Fapar) %>%
    #   mutate(timestamp = timestamp %>% lubridate::ymd() %>% 
    #            as.numeric()) %>% select(timestamp)
    # model_fapar <- splinefun(x=x$timestamp,y=fapar_filled$Fapar,method ='natural')
    model_fapar <- qgam::qgam(Fapar ~ s(timestamp, bs = 'tp',k = k_), data = fapar_filled %>%
                                dplyr::select(timestamp,Fapar) %>%
                                mutate(timestamp = timestamp %>% lubridate::ymd() %>%
                                         as.numeric()), qu = 0.8)
    
    model_lai <- qgam::qgam(Lai ~ s(timestamp, bs = 'tp', k = k_), data = fapar_filled %>%
                              dplyr::select(timestamp,Lai) %>%
                              mutate(timestamp = timestamp %>% lubridate::ymd() %>%
                                       as.numeric()), qu = 0.8)
    # x <- fapar_filled$timestamp %>%lubridate::ymd() %>% as.numeric()
    # y_fapar <- fapar_filled$Fapar
    # y_lai <- fapar_filled$Lai
    # model_fapar <- GauPro::GauPro(x, y_fapar, parallel=TRUE)
    # model_lai <- GauPro::GauPro(x, y_lai, parallel=TRUE)
    # # 
    # pred_FAPAR <- model_fapar$pred(sfn[,"timestamp"] %>%
    #                                  lubridate::as_date() %>%
    #                                  as.numeric())
    # pred_LAI <- model_lai$pred(sfn[,"timestamp"] %>%
    #                              lubridate::as_date() %>%
    #                              as.numeric())
    pred_FAPAR <- predict(model_fapar, tibble(timestamp = sfn[,"timestamp"] %>%
                                                lubridate::as_date() %>%
                                                as.numeric()))
    pred_LAI <- predict(model_lai, tibble(timestamp = sfn[,"timestamp"] %>%
                                            lubridate::as_date() %>%
                                            as.numeric()))
    
    sfn %>% 
      left_join(tibble(timestamp = sfn$timestamp %>% 
                         unique(), 
                       FAPAR = pred_FAPAR) %>% 
                  filter(timestamp > min_data) %>%
                  mutate(FAPAR = case_when(FAPAR<0~0,
                                           FAPAR>1~1,
                                           TRUE~FAPAR)),
                by = "timestamp") %>% 
      left_join(tibble(timestamp = sfn$timestamp %>% 
                         unique(),  LAI = pred_LAI) %>%
                  filter(timestamp > min_data) %>% 
                  mutate(LAI = case_when(LAI<0~0,
                                         LAI>=0~LAI)),
                by = "timestamp") %>% 
      left_join(fapar_filled %>% mutate(timestamp = lubridate::ymd(timestamp)), by = "timestamp")
    
  }else{
    sfn %>%
      left_join(fapar_filled %>% mutate(timestamp = lubridate::ymd(timestamp)), by = "timestamp") %>% 
      mutate(FAPAR = mean(Fapar,na.rm = TRUE),
             LAI = mean(Lai,na.rm = TRUE))
  }
}


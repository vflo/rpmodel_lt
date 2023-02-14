
include_albedo <- function(sfn, albedo){
  albedo_filled <- albedo %>% 
    filter(si_code == unique(sfn$si_code)
    )%>%
    group_by(date,si_code) %>% 
    summarise_all(mean, na.rm = TRUE) %>% 
    dplyr::select(timestamp = date, albedo = Albedo_WSA_vis, albedo_nir = Albedo_WSA_nir) %>%
    mutate(albedo = albedo*0.01,  albedo_nir =  albedo_nir*0.001, source = "MODIS")
  
  min_data = min(albedo_filled$timestamp)
  
  # Build FAPAR model
  k_ <- as.integer(nrow(albedo_filled)*0.05)
  if(k_ >= 3){
    model_albedo <- qgam::qgam(albedo ~ s(timestamp, bs = 'tp',k = k_), data = albedo_filled %>%
                                dplyr::select(timestamp, albedo) %>%
                                mutate(timestamp = timestamp %>% lubridate::ymd() %>%
                                         as.numeric()), qu = 0.5)

    pred_albedo <- predict(model_albedo, tibble(timestamp = sfn[,"timestamp"] %>%
                                                lubridate::as_date() %>%
                                                as.numeric()))
    
    # model_albedo_nir <- qgam::qgam(albedo_nir ~ s(timestamp, bs = 'tp',k = k_), data = albedo_filled %>%
    #                              dplyr::select(timestamp, albedo_nir) %>%
    #                              mutate(timestamp = timestamp %>% lubridate::ymd() %>%
    #                                       as.numeric()), qu = 0.5)
    # 
    # pred_albedo_nir <- predict(model_albedo_nir, tibble(timestamp = sfn[,"timestamp"] %>%
    #                                               lubridate::as_date() %>%
    #                                               as.numeric()))

    
    sfn %>% 
      left_join(tibble(timestamp = sfn$timestamp %>% 
                         unique(), 
                       ALBEDO = pred_albedo) %>% 
                  filter(timestamp > min_data) %>%
                  mutate(ALBEDO = case_when(ALBEDO<0~0,
                                            ALBEDO>1~1,
                                            TRUE~ALBEDO),
                         ),
                by = "timestamp") %>% 
      # left_join(tibble(timestamp = sfn$timestamp %>% 
      #                    unique(), 
      #                  ALBEDO_NIR = pred_albedo_nir) %>% 
      #             filter(timestamp > min_data) %>%
      #             mutate(ALBEDO_NIR = case_when(ALBEDO_NIR<0~0,
      #                                           ALBEDO_NIR>1~1,
      #                                           TRUE~ALBEDO_NIR),
      #             ),
      #           by = "timestamp") %>% 
      left_join(albedo_filled %>% mutate(timestamp = lubridate::ymd(timestamp)), by = "timestamp")
    
  }else{
    sfn %>%
      left_join(albedo_filled %>% mutate(timestamp = lubridate::ymd(timestamp)), by = "timestamp") %>% 
      mutate(ALBEDO = mean(albedo,na.rm = TRUE)#,
             # ALBEDO_NIR = mean(albedo_nir,na.rm = TRUE)
             )
  }
}


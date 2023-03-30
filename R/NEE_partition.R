library(REddyProc)
library(lubridate)
library(tidyverse)
source("R/function_to_get_the_gpp_partition_using_NT_method.R")
set.seed(0815)      # for reproducible results
kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
#### read the files' paths ####
filenames.fluxnet<- list.files(path="R/data/original_sites", "*.csv$", full.names=TRUE,recursive = TRUE)
filename <- filenames.fluxnet[13]
fluxnet_data<-data.table::fread(filename ,  header=T, quote="\"", sep=",",na.strings = "-9999",integer64="character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
site
df <- tibble(timestamp=ind) %>% cbind(as_tibble(fluxnet_data)) #%>% cbind(si_code = site)
plot(df$timestamp,df$USTAR_1_1_1)
#### read metadata ####
sites_metadata <- read_delim("R/data/sites_metadata.csv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
#Sites where the nee partition failed:
#"CR-SoC" "PR-xGU" "US-xSP" "US-xTE"
#### Estimate partition ####
purrr::map2(.x=as.list(filenames.fluxnet)[3],
            .y=sites_metadata[3,]%>%
              split(seq(nrow(.))),
            .f=function(x = .x, y=.y){
              site <- do.call(rbind, strsplit(basename(x), "_"))[,2]
              if(site != y$site){stop()}
              part_nee(x,y$lat, y$long,"R/data/sites",
                       start_q =NA,end_q = NA,dts=48,
                       offset=y$UTC_OFFSET)
              }
            )







#### Obtain final database ####
filenames_fluxnet <- list.files(path="R/data/sites", "*.csv$", full.names=TRUE,recursive = TRUE)
filenames_fluxnet_name <- list.files(path="R/data/sites", "*.csv$", full.names=FALSE,recursive = TRUE)
index <- 12
purrr::map(as.list(c(10,12,13,29)),function(index){
  filename <- filenames_fluxnet[index]
  filename_name <- filenames_fluxnet_name[index]
  fluxnet_data<-data.table::fread(filename ,  header=T, quote="\"", sep=",")
  ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M",
                  tz = "GMT") %>% as.POSIXct()
  time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
  t_conv_f <- 3600 * time.freq
  site <- do.call(rbind, strsplit(basename(filename), "_"))[,4]
  site
  df <- tibble(DateTime=ind) %>% cbind(as_tibble(fluxnet_data)) #%>% cbind(si_code = site)
  names(df)
  
  df <- df %>% dplyr::select(-starts_with("WS_MAX"),
                             -starts_with("CO2_M"),
                             -starts_with("H2O"),
                             -starts_with("H_S")) %>% 
    mutate(T_CANOPY = ifelse(length(grep("^T_CANOPY$", names(df), value = TRUE))==0,NA_real_,
                             get(grep("^T_CANOPY", names(df), value = TRUE)[1]))) %>% #Some sites don't have T_CANOPY variable
    mutate(DateTime = DateTime,
           TIMESTAMP_START = TIMESTAMP_START,
           TIMESTAMP_END = TIMESTAMP_END,
           Tair = TA,
           Ustar = USTAR,
           rH = RH,
           NEE = NEE_f,
           GPP = GPP,#
           RECO = RECO,
           H = case_when(length(grep("^H", names(df), value = TRUE))>0~
                           rowMeans(dplyr::select(df,starts_with("H")),na.rm=TRUE),
                         TRUE~NA_real_),
           LW_OUT = case_when(length(grep("^LW_OUT", names(df), value = TRUE))>0~
                           rowMeans(dplyr::select(df,starts_with("LW_OUT")),na.rm=TRUE),
                         TRUE~NA_real_),
           SW_IN = SW_IN,
           VPD = VPD_PI,
           PPFD = SW_IN * (kfFEC*(1 - kalb_vis)),#bigleaf::Rg.to.PPFD(SW_IN, J_to_mol = 4.6, frac_PAR = 0.5),#rowMeans(select(df,starts_with("PPFD_IN")),na.rm=TRUE),
           WS = rowMeans(dplyr::select(df,starts_with("WS")),na.rm=TRUE),
           NETRAD = NETRAD,
           CO2 = case_when(length(grep("^CO2_PI$", names(df), value = TRUE))>0~
                             rowMeans(dplyr::select(df,starts_with("CO2_PI")),na.rm=TRUE),
                           TRUE~rowMeans(dplyr::select(df,-starts_with("CO2_PI")) %>% 
                                           dplyr::select(starts_with("CO2")),na.rm=TRUE)),
           PA = case_when(length(grep("^PA_PI$", names(df), value = TRUE))>0~
                            rowMeans(dplyr::select(df,starts_with("PA_PI")),na.rm=TRUE),
                          TRUE~rowMeans(dplyr::select(df,-starts_with("PA_PI")) %>% 
                                          dplyr::select(starts_with("PA")),na.rm=TRUE)),
           LE = case_when(length(grep("^LE_PI$", names(df), value = TRUE))>0~
                            rowMeans(dplyr::select(df,starts_with("LE_PI")),na.rm=TRUE),
                          TRUE~rowMeans(dplyr::select(df,-starts_with("LE_PI")) %>% 
                                          dplyr::select(starts_with("LE")),na.rm=TRUE)),
           P = case_when(length(grep("^P_PI$", names(df), value = TRUE))>0~
                           rowMeans(dplyr::select(df,starts_with("P_PI")),na.rm=TRUE),
                         TRUE~rowMeans(dplyr::select(df,-starts_with(c("PP","P_PI","PA"))) %>% 
                                         dplyr::select(starts_with("P")),na.rm=TRUE)),
           SWC = case_when(length(grep("^SWC_PI$", names(df), value = TRUE))>0~
                             rowMeans(dplyr::select(df,starts_with("SWC_PI")),na.rm=TRUE),
                           TRUE~rowMeans(dplyr::select(df,-starts_with("SWC_PI")) %>% 
                                           dplyr::select(starts_with("SWC")),na.rm=TRUE)),
           # Tcan =  rowMeans(select(df,starts_with("T_CANOPY")),na.rm=TRUE)
           Tcan =  case_when(length(grep("^T_CANOPY", names(.), value = TRUE))>0~
                               get(grep("^T_CANOPY", names(.), value = TRUE)[1]),
                             TRUE~NA_real_)
           ) %>%
             dplyr::select(DateTime,TIMESTAMP_START, TIMESTAMP_END, Tair, Ustar, rH,NEE, GPP,RECO,H,
                           LW_OUT, SW_IN, VPD, PPFD, WS, NETRAD, CO2, PA, LE, P, SWC, Tcan)
  
  summary(df)
  
  write.csv(df,paste0("R/data/final_sites/",filename_name), row.names=FALSE)

})




#### Include canopy temperatures to HA1, Me2, NR1, Wrc ####
filenames_fluxnet <- list.files(path="R/data/final_sites", "*.csv$", full.names=TRUE,recursive = TRUE)
filenames_fluxnet

#Ha1
index <- 10
filename <- filenames_fluxnet[index]
fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",")
load(file="R/data/HF.RData")
names(HF)
foo <- HF %>% 
  mutate(Tcan = rowMeans(dplyr::select(.,ACRU_mean, BEPA_mean, PIST_mean, QURU_mean), na.rm = TRUE),
         DateTime = (date - lubridate::hours(8)-minutes(30) ) %>% force_tz(tzone = "UTC")
         ) %>% 
  dplyr::select(DateTime,Tcan,ACRU_mean, BEPA_mean, PIST_mean, QURU_mean)
fluxnet_data <- fluxnet_data %>% 
  dplyr::select(-Tcan) %>% 
  left_join(foo )
write.csv(fluxnet_data, filename, row.names=FALSE)

#Me2
index <- 12
filename <- filenames_fluxnet[index]
fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",")
load(file="R/data/MR.RData")
names(MR)
foo <- MR %>% 
  mutate(Tcan = Tcan_Avg_corr,2,
         DateTime = (date - lubridate::minutes(15))%>% force_tz(tzone = "UTC")) %>% 
  dplyr::select(DateTime,Tcan)
fluxnet_data <- fluxnet_data %>% 
  mutate(DateTime = DateTime) %>% 
  dplyr::select(-Tcan) %>% 
  left_join(foo)
write.csv(fluxnet_data, filename, row.names=FALSE)

#NR1
index <- 13
filename <- filenames_fluxnet[index]
fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",")
load(file="R/data/NW.Rdata")
names(NW)
foo <- NW %>% 
  mutate(Tcan = Tcan,
         DateTime = (lubridate::ymd_hms(date)+ lubridate::hours(9) - lubridate::minutes(15)) %>% force_tz(tzone = "UTC")) %>% 
  dplyr::select(DateTime,Tcan)
fluxnet_data <- fluxnet_data %>% 
  dplyr::select(-Tcan) %>% 
  left_join(foo)
write.csv(fluxnet_data, filename, row.names=FALSE)

#WRC
index <- 29
filename <- filenames_fluxnet[index]
fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",")
load(file="R/data/WR.RData")
names(WR)
foo <- WR %>% 
  mutate(Tcan = Tcan_Avg_corr,
         DateTime = (lubridate::ymd_hms(date)- lubridate::hours(8)) %>% force_tz(tzone = "UTC")) %>% 
  dplyr::select(DateTime,Tcan)
fluxnet_data <- fluxnet_data %>% 
  dplyr::select(-Tcan) %>% 
  left_join(foo)
write.csv(fluxnet_data, filename, row.names=FALSE)


# filenames.fluxnet
# filename <- filenames.fluxnet[7]
# fluxnet_data <- data.table::fread(filename, header = T, 
#                                   quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
#                                   colClasses = "numeric", integer64 = "character")
# ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M", 
#                 tz = "GMT") %>% as.POSIXct()
# time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
# t_conv_f <- 3600 * time.freq
# site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
# site
# df <- tibble(timestamp=ind) %>% cbind(as_tibble(fluxnet_data)) #%>% cbind(si_code = site)
# 
# names(df)
# df <- df %>%
#   mutate(DateTime = timestamp,
#          TIMESTAMP_START = TIMESTAMP_START,
#          TIMESTAMP_END = TIMESTAMP_END,
#          Tair = rowMeans(select(df,starts_with("TA_")),na.rm=TRUE),
#          Ustar = USTAR,
#          rH = rowMeans(select(df,starts_with("RH")),na.rm=TRUE),
#          NEE = NEE_PI,
#          GPP = NA,#GPP_PI,
#          FC = rowMeans(select(df,starts_with("FC")),na.rm=TRUE),
#          H = NA,#H_2_1_1,
#          Rg = rowMeans(select(df,starts_with("SW_IN")),na.rm=TRUE),
#          VPD = rowMeans(select(df,starts_with("VPD")),na.rm=TRUE),
#          PPFD = bigleaf::Rg.to.PPFD(Rg, J_to_mol = 4.6, frac_PAR = 0.5),#rowMeans(select(df,starts_with("PPFD_IN")),na.rm=TRUE),
#          LEAF_WET = NA,#rowMeans(select(df,starts_with("LEAF_WET")),na.rm=TRUE),
#          WS = WS_1_1_1,#rowMeans(select(df,starts_with(paste0("WS_",c(1,2)))),na.rm=TRUE),
#          Rnet = rowMeans(select(df,starts_with("NETRAD")),na.rm=TRUE),
#          CO2 = rowMeans(select(df,starts_with("CO2")),na.rm=TRUE),
#          PA = rowMeans(select(df,starts_with("PA")),na.rm=TRUE),
#          LE = rowMeans(select(df,starts_with("LE")),na.rm=TRUE),
#          P = rowMeans(select(df,starts_with("P")),na.rm=TRUE),
#          SWC = rowMeans(select(df,starts_with("SWC_")),na.rm=TRUE),
#          Tcan =  rowMeans(select(df,starts_with("T_CANOPY")),na.rm=TRUE)) %>%
#            dplyr::select(DateTime,TIMESTAMP_START, TIMESTAMP_END, Tair, Ustar, rH,NEE,
#                          Rg, VPD, PPFD, LEAF_WET, WS, Rnet, CO2, PA, LE, P, SWC, Tcan,GPP,FC,H)
# 
# #Only one variable
# # df <- df %>%
# #   mutate(DateTime = timestamp,
# #          Tair = TA,
# #          Ustar = USTAR,
# #          rH = RH,
# #          NEE = NEE_PI,
# #          Rg = SW_IN,
# #          VPD = VPD_PI,
# #          PPFD = bigleaf::Rg.to.PPFD(Rg, J_to_mol = 4.6, frac_PAR = 0.5),
# #          LEAF_WET = NA,
# #          WS = WS,
# #          Rnet = NETRAD,
# #          CO2 = CO2,
# #          PA = PA,
# #          LE = LE,
# #          P = P,
# #          SWC = SWC,
# #          Tcan =  NA) %>%
# #   dplyr::select(DateTime, Tair, Ustar, rH,NEE, Rg, VPD, PPFD, LEAF_WET, WS, Rnet, CO2, PA, LE, P, SWC, Tcan) 
# 
# #NEE not present in the data
# names(df)
# df <- df %>%
#   mutate(DateTime = timestamp,
#          TIMESTAMP_START = TIMESTAMP_START,
#          TIMESTAMP_END = TIMESTAMP_END,
#          TA = rowMeans(select(df,starts_with("TA_")),na.rm=TRUE),
#          USTAR = USTAR_2_1_1,
#          RH = rowMeans(select(df,starts_with("RH_")),na.rm=TRUE),
#          NEE = NA,#NEE_PI,
#          GPP = NA,#GPP_PI,
#          FC = rowMeans(select(df,starts_with("FC")),na.rm=TRUE),
#          H = H_2_1_1,
#          SW_IN = rowMeans(select(df,starts_with("SW_IN")),na.rm=TRUE),
#          VPD = rowMeans(select(df,starts_with("VPD")),na.rm=TRUE),
#          PPFD = bigleaf::Rg.to.PPFD(SW_IN, J_to_mol = 4.6, frac_PAR = 0.5),#rowMeans(select(df,starts_with("PPFD_IN")),na.rm=TRUE),
#          LEAF_WET = NA,#rowMeans(select(df,starts_with("LEAF_WET")),na.rm=TRUE),
#          WS = WS_1_1_1,#rowMeans(select(df,starts_with(paste0("WS_",c(1,2)))),na.rm=TRUE),
#          Rnet = rowMeans(select(df,starts_with("NETRAD")),na.rm=TRUE),
#          CO2 = rowMeans(select(df,starts_with("CO2")),na.rm=TRUE),
#          PA = rowMeans(select(df,starts_with("PA")),na.rm=TRUE),
#          LE = rowMeans(select(df,starts_with("LE")),na.rm=TRUE),
#          P = rowMeans(select(df,starts_with("P")),na.rm=TRUE),
#          SWC = rowMeans(select(df,starts_with("SWC_")),na.rm=TRUE),
#          Tcan =  rowMeans(select(df,starts_with("T_CANOPY")),na.rm=TRUE)) %>%
#   dplyr::select(DateTime,TIMESTAMP_START, TIMESTAMP_END, TA, USTAR, RH, NEE,
#                 SW_IN, VPD, PPFD, LEAF_WET, WS, Rnet, CO2, PA, LE, P, SWC, Tcan,GPP,FC,H)
# foo <- REddyProc::read_from_ameriflux22(df)
# 
# names(foo)
# df <- foo %>% cbind(df %>% dplyr::select(PPFD,LEAF_WET,WS,Rnet,CO2,PA,P,SWC,Tcan))
# df <- df %>%
#   slice(56641:76800) %>%
#   rowwise() %>%
#   mutate(Rg = case_when(Rg<0~0,TRUE~Rg),
#          NEE = case_when(as.numeric(NEE)<=100~NEE,
#                          TRUE~NaN))
# # 
# # 
# # 
# # df <- filterLongRuns(df, "NEE")
# # 
# # 
# # summary(df)
# # 
# # 
# # df$VPD <- fCalcVPDfromRHandTair(df$rH, df$Tair)
# EProc <- sEddyProc$new(site, df, c('NEE','Rg','Tair','VPD', 'Ustar'))
# EProc$sSetLocationInfo(LatDeg = sites_metadata$lat[3], LongDeg = sites_metadata$long[3], TimeZoneHour = -3)  #Location
# # 	
# # 
# # seasonStarts <- as.data.frame( do.call( rbind, list(
# #   c(100,2014)
# #   # ,c(70,2017)
# #   # ,c(200,2017)
# #   # ,c(320,2017)
# #   ,c(200,2018)
# #   # ,c(305,2018)
# #   # ,c(200,2019)
# #   # ,c(305,2019)
# #   # ,c(200,2020)
# #   # ,c(305,2020)
# #   # ,c(200,2021)
# #   # ,c(305,2021)
# # )))
# # seasonFactor <- usCreateSeasonFactorYdayYear(
# #   df$DateTime - 15*60, starts = seasonStarts)
# # plot( NEE ~ DateTime, df)
# # plot( Tcan ~ DateTime, df)
# # plot( Ustar ~ DateTime, df)
# # seasonStartsDate <- fConvertTimeToPosix( data.frame(Year = seasonStarts[,2]
# #                                                     , DoY = seasonStarts[,1], Hour = 0.25), 'YDH'
# #                                          , Year = "Year", Day = "DoY", Hour = "Hour")
# # abline( v = seasonStartsDate$DateTime)
# # 
# EProc$sEstimateUstarScenarios( nSample = 100L, probs = c(0.05,0.5,0.95))#,seasonFactor = seasonFactor)
# EProc$useSeaonsalUStarThresholds()
# EProc$sGetUstarScenarios()
# # 
# EProc$sMDSGapFillUStarScens('NEE', FillAll = TRUE)
# EProc$sMDSGapFill('NEE', FillAll.b=TRUE)
# # 
# grep("^NEE.*_f$", colnames( EProc$sExportResults()), value = TRUE )
# EProc$sPlotFingerprintY('NEE_f', Year = 2018)
# 
# EProc$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)
# EProc$sMDSGapFill('VPD', FillAll = FALSE,  minNWarnRunLength = NA)
# EProc$sMDSGapFill('Rg', FillAll = FALSE,  minNWarnRunLength = NA)
# EProc$sFillVPDFromDew() # fill longer gaps still present in VPD_f
# EProc$sMRFluxPartitionUStarScens()
# EProc$sGLFluxPartitionUStarScens()
# EProc$sMRFluxPartition(TempVar = "Tair", QFTempVar = "NEE_fqc", 
#                        RadVar = "Rg", parsE0Regression=list(TempRange=1L))
# EProc$sGLFluxPartition(TempVar = "Tair", QFTempVar = "NEE_fqc", 
#                        RadVar = "Rg", VPDVar ='VPD',
#                        controlGLPart = partGLControl(nBootUncertainty=0L, 
#                                                      isAssociateParmsToMeanOfValids=FALSE , isLasslopPriorsApplied=TRUE, isBoundLowerNEEUncertainty=FALSE , smoothTempSensEstimateAcrossTime=FALSE,useNightimeBasalRespiration=T))
# 
# grep("^GPP.*_f$", colnames( EProc$sExportResults()), value = TRUE )
# EProc$sPlotFingerprintY('GPP_DT_U50',#'GPP_U50_f',
#                         Year = 2018)
# 
# REddyProc::sEddyProc_sGLFluxPartition(FilledEddyData)
# 
# FilledEddyData <- EProc$sExportResults()
# uStarSuffixes <- colnames(EProc$sGetUstarScenarios())[-1]
# GPPAggCO2 <- sapply( uStarSuffixes, function(suffix) {
#   GPPHalfHour <- FilledEddyData[[paste0("GPP_",suffix,"_f")]]
#   mean(GPPHalfHour, na.rm = TRUE)
# })
# molarMass <- 12.011
# GPPAgg <- GPPAggCO2 * 1e-6 * molarMass * 3600*24*365.25
# print(GPPAgg)
# 
# (max(GPPAgg) - min(GPPAgg)) / median(GPPAgg)
# 
# 
# FilledEddyData <- EProc$sExportResults()
# CombinedData <- cbind(df, FilledEddyData)
# # CombinedData <- df
# # fWriteDataframeToFile(CombinedData, paste0(site,'-Results.txt'), Dir = "R/data/sites")
# # 
# # 
# # # 
# # # 
# # # 
# # # 
# # # url <- "http://cran.r-project.org/src/contrib/Archive/mlegp/mlegp_3.1.8.tar.gz"
# # # pkgFile <- "mlegp_3.1.8.tar.gz"
# # # download.file(url = url, destfile = pkgFile)
# # # 
# # # # Expand the zip file using whatever system functions are preferred
# # # 
# # # # look at the DESCRIPTION file in the expanded package directory
# # # 
# # # # Install dependencies list in the DESCRIPTION file
# # # 
# # # install.packages(c("ada", "ipred", "evd"))
# # # 
# # # # Install package
# # # install.packages(pkgs=pkgFile, type="source", repos=NULL)

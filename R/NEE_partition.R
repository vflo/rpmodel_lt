library(REddyProc)
library(lubridate)
library(tidyverse)
set.seed(0815)      # for reproducible results

#### read the files' paths
filenames.fluxnet<- list.files(path="R/data", "*.csv$", full.names=TRUE,recursive = TRUE)
filename <- filenames.fluxnet[3]
fluxnet_data <- data.table::fread(filename, header = T, 
                                  quote = "\"", sep = ",", na.strings = c("NA", "-9999"), 
                                  colClasses = "numeric", integer64 = "character")
ind <- strptime(fluxnet_data$TIMESTAMP_START, format = "%Y%m%d%H%M", 
                tz = "GMT") %>% as.POSIXct()
time.freq <- abs(as.numeric(ind[1] - ind[2], units = "hours"))
t_conv_f <- 3600 * time.freq
site <- do.call(rbind, strsplit(basename(filename), "_"))[,2]
df <- tibble(timestamp=ind) %>% cbind(as_tibble(fluxnet_data)) %>% cbind(si_code = site)


df <- df %>%
  mutate(DateTime = timestamp,
         Tair = rowMeans(select(df,starts_with("TA_")),na.rm=TRUE),
         Ustar = USTAR,
         rH = RH,
         NEE = NEE_PI,
         Rg = rowMeans(select(df,starts_with("SW_IN")),na.rm=TRUE),
         VPD = rowMeans(select(df,starts_with("VPD")),na.rm=TRUE),
         PPFD = rowMeans(select(df,starts_with("PPFD_IN")),na.rm=TRUE),
         LEAF_WET = rowMeans(select(df,starts_with("LEAF_WET")),na.rm=TRUE),
         WS = rowMeans(select(df,starts_with(paste0("WS_",c(1,2)))),na.rm=TRUE),
         Rnet = rowMeans(select(df,starts_with("NETRAD")),na.rm=TRUE),
         CO2 = rowMeans(select(df,starts_with("CO2")),na.rm=TRUE),
         PA = rowMeans(select(df,starts_with("PA_")),na.rm=TRUE),
         LE = rowMeans(select(df,starts_with("LE_")),na.rm=TRUE),
         P = P,
         SWC = rowMeans(select(df,starts_with("SWC_")),na.rm=TRUE),
         Tcan =  rowMeans(select(df,starts_with("T_CANOPY")),na.rm=TRUE)) %>%
           dplyr::select(DateTime, Tair, Ustar, rH,NEE, Rg, VPD, PPFD, LEAF_WET, WS, Rnet, CO2, PA, LE, P, SWC, Tcan) 

df <- filterLongRuns(df, "NEE")


summary(df)


df$VPD <- fCalcVPDfromRHandTair(df$rH, df$Tair)
EProc <- sEddyProc$new('US-xDS', df, c('NEE','Rg','Tair','VPD', 'Ustar'))
EProc$sSetLocationInfo(LatDeg = 28.1250, LongDeg = -81.4362, TimeZoneHour = -5)  #Location of Gebesee

EProc$sEstimateUstarScenarios( nSample = 100L, probs = c(0.05,0.5,0.95))
EProc$useSeaonsalUStarThresholds()
EProc$sGetUstarScenarios()

EProc$sMDSGapFillUStarScens('NEE', FillAll = FALSE)

grep("^NEE.*_f$", colnames( EProc$sExportResults()), value = TRUE )
EProc$sPlotFingerprintY('NEE_U50_f', Year = 2020)

EProc$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)     
EProc$sMDSGapFill('VPD', FillAll = FALSE,  minNWarnRunLength = NA)     
EProc$sFillVPDFromDew() # fill longer gaps still present in VPD_f
EProc$sMRFluxPartitionUStarScens()
grep("^GPP.*_f$", colnames( EProc$sExportResults()), value = TRUE )
EProc$sPlotFingerprintY('GPP_U50_f', Year = 2021)


FilledEddyData <- EProc$sExportResults()
uStarSuffixes <- colnames(EProc$sGetUstarScenarios())[-1]
GPPAggCO2 <- sapply( uStarSuffixes, function(suffix) {
  GPPHalfHour <- FilledEddyData[[paste0("GPP_",suffix,"_f")]]
  mean(GPPHalfHour, na.rm = TRUE)
})
molarMass <- 12.011
GPPAgg <- GPPAggCO2 * 1e-6 * molarMass * 3600*24*365.25
print(GPPAgg)

(max(GPPAgg) - min(GPPAgg)) / median(GPPAgg) 


FilledEddyData <- EProc$sExportResults()
CombinedData <- cbind(df, FilledEddyData)
fWriteDataframeToFile(CombinedData, 'US-xDS-Results.txt', Dir = "R/data/sites")

#function to get the gpp partition using NT method
part_nee<-function(filename,lat,lon,outpath=getwd(),start_q = NA, end_q = NA, dts=48, offset = NA){
	#filename="X:/home/WORK/data_input/europaflux/SE-Svb_2014-2016_L2_.csv";lat=avail.stations$lat[63];lon= avail.stations$lon[63]
	fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",",na.strings = "-9999",integer64="character")
	#eur
	# site<-do.call(rbind,strsplit(basename(filename),'_'))[,2]
	#ameri
	site<-do.call(rbind,strsplit(basename(filename),'_'))[,2]
	kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
	kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
	#############################################################################################################
	###1. organize the df
	#############################################################################################################
	#choose NEE column if null, choose FC (carbon flux column)
	if(length(grep("NEE", names(fluxnet_data), value = TRUE))==0 & length(grep("FC", names(fluxnet_data), value = TRUE))>=1){
		namec<-grep("FC", names(fluxnet_data), value = TRUE)[1]
		fluxnet_data$NEE_f<-fluxnet_data[,..namec]
		fluxnet_data$NEE_fqc=rep(0,length(fluxnet_data[,1])) 
	}else{
		namec<-grep("NEE", names(fluxnet_data), value = TRUE)[1]
		fluxnet_data$NEE_f<-fluxnet_data[,..namec]
		fluxnet_data$NEE_fqc=rep(0,length(fluxnet_data[,1])) 
	}
	#some NEE_PI report all NAs..??
	if(sum(is.na(fluxnet_data$NEE_f))==length(fluxnet_data$NEE_f)){
		namec<-grep("FC", names(fluxnet_data), value = TRUE)[1]
		fluxnet_data$NEE_f<-fluxnet_data[,..namec]
	}
	
	#some NEE_PI report some NAs..??
	if(length(grep("FC", names(fluxnet_data), value = TRUE))>=1){
	  fluxnet_data <- fluxnet_data %>% 
	    mutate(NEE_f = coalesce(NEE_f,get(grep("FC", names(fluxnet_data), value = TRUE)[1])))
	}
	
	#add storage flux
	if(length(grep("SC", names(fluxnet_data), value = TRUE))>=1){
	  namec<-grep("SC", names(fluxnet_data), value = TRUE)[1]
	  fluxnet_data$NEE_f<-fluxnet_data$NEE_f+fluxnet_data[,..namec]
	}
	
	
	##### correct the names
	if(length(grep("^PPFD_IN_PI", names(fluxnet_data), value = TRUE))>0){
	  nameppfd<-grep("^PPFD_IN_PI", names(fluxnet_data), value = TRUE)[1]
	  fluxnet_data$PPFD_IN<-rowMeans(fluxnet_data[,..nameppfd],na.rm=T)
	}
	
	#some stattions report ony PPFD and not SW_IN
	if((length(grep("SW_IN", names(fluxnet_data), value = TRUE))==0 & length(grep("PPFD_IN", names(fluxnet_data), value = TRUE))>=1)){
	  nameppfd<-grep("PPFD_IN", names(fluxnet_data), value = TRUE)
	  if(any(nameppfd == "PPFD_IN")){nameppfd <- "PPFD_IN"}else{nameppfd<-nameppfd[1]}
	  fluxnet_data$SW_IN<-fluxnet_data[,..nameppfd]/(kfFEC*(1 - kalb_vis))
	  fluxnet_data$SW_IN[fluxnet_data$SW_IN<0]<-0
	}else if(length(grep("SW_IN", names(fluxnet_data), value = TRUE))>=1 & sum(is.na(fluxnet_data$SW_IN))==length(fluxnet_data$SW_IN)&length(grep("PPFD_IN", names(fluxnet_data), value = TRUE))>=1){
	  nameppfd<-grep("PPFD_IN", names(fluxnet_data), value = TRUE)
	  if(any(nameppfd == "PPFD_IN")){nameppfd <- "PPFD_IN"}else{nameppfd<-nameppfd[1]}
	  fluxnet_data$SW_IN<-fluxnet_data[,..nameppfd]/(kfFEC*(1 - kalb_vis))
	  fluxnet_data$SW_IN[fluxnet_data$SW_IN<0]<-0
	}
	
	##sw_in
	if(length(grep("^SW_IN$", names(fluxnet_data), value = TRUE))==0){
		namesw<-grep("^SW_IN", names(fluxnet_data), value = TRUE)[1]
		fluxnet_data$SW_IN<-fluxnet_data[,..namesw]
		
	}
	

	
	##USTAR
	if(length(grep("^USTAR$", names(fluxnet_data), value = TRUE))==0){
		nameneustr<-grep("^USTAR", names(fluxnet_data), value = TRUE)[1]
		fluxnet_data$USTAR<-fluxnet_data[,..nameneustr]
		
	}
	# ##VPD
	# if(length(grep("^VPD$", names(fluxnet_data), value = TRUE))==0){
	# 	namenevpd<-grep("^VPD", names(fluxnet_data), value = TRUE)[1]
	# 	fluxnet_data$VPD_PI<-fluxnet_data[,..namenevpd]
	# 	
	# }

	
	##netrad
	if(length(grep("^NETRAD$", names(fluxnet_data), value = TRUE))>0){
  	if(all(is.na(fluxnet_data$NETRAD))&!all(is.na(fluxnet_data$SW_IN),is.na(fluxnet_data$SW_OUT),
  	                                        is.na(fluxnet_data$LW_IN),is.na(fluxnet_data$LW_OUT))){
  	  fluxnet_data$NETRAD<- fluxnet_data$SW_IN - fluxnet_data$SW_OUT + fluxnet_data$LW_IN - fluxnet_data$LW_OUT
  	}
	 }
	if(length(grep("^NETRAD$", names(fluxnet_data), value = TRUE))==0){
	  namenetr<-grep("^NETRAD", names(fluxnet_data), value = TRUE)[1]
	  fluxnet_data$NETRAD<-fluxnet_data[,..namenetr]
	  if(all(is.na(fluxnet_data$NETRAD))&!all(is.na(fluxnet_data$SW_IN),is.na(fluxnet_data$SW_OUT),
	                                          is.na(fluxnet_data$LW_IN),is.na(fluxnet_data$LW_OUT))){
	    fluxnet_data$NETRAD<- fluxnet_data$SW_IN - fluxnet_data$SW_OUT + fluxnet_data$LW_IN - fluxnet_data$LW_OUT
	  }
	}
	
	####get TA average
	if(length(grep("^TA$", names(fluxnet_data), value = TRUE))==0){
		nameta<-grep("^TA_", names(fluxnet_data), value = TRUE)[1]
		# fluxnet_data$TA<-rowMeans(fluxnet_data[,..nameta],na.rm=T)
		fluxnet_data$TA<-fluxnet_data[,..nameta]
	}
	
	if(length(grep("^TA_PI", names(fluxnet_data), value = TRUE))>0){
	  nameta<-grep("^TA_PI", names(fluxnet_data), value = TRUE)[1]
	  # fluxnet_data$TA<-rowMeans(fluxnet_data[,..nameta],na.rm=T)
	  fluxnet_data$TA<-fluxnet_data[,..nameta]
	}
	
	####get RH average
	if(length(grep("^RH$", names(fluxnet_data), value = TRUE))==0){
		namerh<-grep("^RH_", names(fluxnet_data), value = TRUE)[1]
		# fluxnet_data$RH<-rowMeans(fluxnet_data[,..namerh],na.rm=T)
		fluxnet_data$RH<-fluxnet_data[,..namerh]
		fluxnet_data$RH[fluxnet_data$RH>100]<-99.99
		
	}
	#calc vpd hPa if not provided
	if(length(grep("^VPD", names(fluxnet_data), value = TRUE))==0){
		fluxnet_data$VPD_PI <- fCalcVPDfromRHandTair(fluxnet_data$RH, fluxnet_data$TA)
		
	}else if(length(grep("^VPD", names(fluxnet_data), value = TRUE))>=1){
		namevpd<-grep("^VPD_", names(fluxnet_data), value = TRUE)[1]
		# fluxnet_data$VPD_PI<-rowMeans(fluxnet_data[,..namevpd],na.rm=T)
		fluxnet_data$VPD_PI<-fluxnet_data[,..namevpd]
		#fluxnet_data$VPD_PI<-fluxnet_data$VPD
	}
	fluxnet_data$VPD_PI[fluxnet_data$VPD_PI<0]<-0.0001

	#############################################################################################################
	###2. build the object required by REddyProc
	#############################################################################################################	
	
  fluxnet_data$Tair <- fluxnet_data$TA
	fluxnet_data$Ustar <- fluxnet_data$USTAR
	fluxnet_data$Rnet <- fluxnet_data$NETRAD
	fluxnet_data$Rg <- fluxnet_data$SW_IN
	fluxnet_data$Rg[fluxnet_data$Rg<0] <- 0
	fluxnet_data$VPD <- fluxnet_data$VPD_PI
	fluxnet_data$NEE <- fluxnet_data$NEE_f
	
	ind<-as.POSIXct(strptime(fluxnet_data$TIMESTAMP_START,format="%Y%m%d%H%M",tz="GMT"))
	yrs<-unique(format(ind,'%Y'))
	fluxnet_data$DateTime<-ind
	if(all(!is.na(start_q),!is.na(start_q))){
	  fluxnet_data <- fluxnet_data %>% filter(TIMESTAMP_START>start_q, TIMESTAMP_END<=end_q)
	}
	# fluxnet_data[which(fluxnet_data$NEE_f > 100),'NEE_f'] <- NA

	EddyProc.C <- sEddyProc$new(site,fluxnet_data,DTS=dts)
	## gap fill NEE
	EddyProc.C$sMDSGapFill('NEE', FillAll=TRUE)
	## get time zone to calc potential rad thus NT
	if(is.na(offset)){
  	if(lon<0){
  		offset = floor(lon * 24 / 360)
  	}else{
  		offset =  floor(lon * 24 / 360)
  	}
	}
	 

	EddyProc.C$sSetLocationInfo(LatDeg=lat, LongDeg=lon, TimeZoneHour=offset)
	# EddyProc.C$sEstimateUstarScenarios( nSample = 100L, probs = c(0.05,0.5,0.95),
	#                                     ctrlUstarEst = usControlUstarEst(minValidBootProp = 0.3),
	#                                     ctrlUstarSub = usControlUstarSubsetting(taClasses = 4, 
	#                                                                             UstarClasses = 10, 
	#                                                                             swThr = 10, 
	#                                                                             minRecordsWithinTemp = 50, 
	#                                                                             minRecordsWithinSeason = 80, 
	#                                                                             minRecordsWithinYear = 1500)
	#                                     )#,seasonFactor = seasonFactor)
	# 
	# 
	# EddyProc.C$useSeaonsalUStarThresholds()
	# # foo <- EddyProc.C$sGetUstarScenarios()
	# EddyProc.C$sGetUstarScenarios()
	# EddyProc.C$sMDSGapFillUStarScens('NEE', FillAll = FALSE)
	# EddyProc.C$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)
	# EddyProc.C$sMDSGapFill('VPD', FillAll = FALSE,  minNWarnRunLength = NA)
	# EddyProc.C$sMDSGapFill('Rg', FillAll = FALSE,  minNWarnRunLength = NA)
	# EddyProc.C$sFillVPDFromDew() # fill longer gaps still present in VPD_f
	#get GPP NT
	# EddyProc.C$sMRFluxPartitionUStarScens( parsE0Regression=list(TempRange=3L))
	EddyProc.C$sMRFluxPartition(TempVar = "Tair", QFTempVar = "NEE_fqc", RadVar = "Rg",
	                            parsE0Regression=list(TempRange=1L))
	
	#if GPP NT did not work, get GPP_DT
	if(is.null(EddyProc.C@.xData$sTEMP$GPP_f)){
		EddyProc.C@.xData$sDATA$Rg_fqc<-EddyProc.C@.xData$sDATA$NEE_fqc
		EddyProc.C@.xData$sDATA$NEE_fsd<-EddyProc.C@.xData$sTEMP$NEE_f_fsd
		EddyProc.C$sGLFluxPartition(TempVar = "Tair", QFTempVar = "NEE_fqc", RadVar = "Rg", controlGLPart = partGLControl(nBootUncertainty=0L, isAssociateParmsToMeanOfValids=FALSE , isLasslopPriorsApplied=TRUE, isBoundLowerNEEUncertainty=FALSE , smoothTempSensEstimateAcrossTime=FALSE,useNightimeBasalRespiration=T))
		# EddyProc.C$sGLFluxPartitionUStarScens()
		fluxnet_data$GPP<-EddyProc.C@.xData$sTEMP$GPP_DT
		fluxnet_data$RECO<-EddyProc.C@.xData$sTEMP$Reco_DT
		# fluxnet_data$GPP_05 <- EddyProc.C@.xData$sTEMP$GPP_DT_U05_f
		# fluxnet_data$GPP_95 <- EddyProc.C@.xData$sTEMP$GPP_DT_U95_f
		# fluxnet_data$SW_IN_POT<-EddyProc.C@.xData$sTEMP$PotRad
		fluxnet_data$DateTime <- NULL
		#write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_DT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
		write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_DT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
		
	}else{
		fluxnet_data$GPP<-EddyProc.C@.xData$sTEMP$GPP_f
		fluxnet_data$RECO<-EddyProc.C@.xData$sTEMP$Reco
		# fluxnet_data$GPP_05 <- EddyProc.C@.xData$sTEMP$GPP_U05_f
		# fluxnet_data$GPP_95 <- EddyProc.C@.xData$sTEMP$GPP_U95_f
		# fluxnet_data$SW_IN_POT<-EddyProc.C@.xData$sTEMP$PotRad_NEW
		fluxnet_data$DateTime <- NULL
		#write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_NT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
		write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_NT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
	}
	
	
}


# 
# #function to get the gpp partition using NT method
# part_nee_alt<-function(filename,lat,lon,outpath=getwd(),start_q = NA, end_q = NA){
#   #filename="X:/home/WORK/data_input/europaflux/SE-Svb_2014-2016_L2_.csv";lat=avail.stations$lat[63];lon= avail.stations$lon[63]
#   fluxnet_data<-data.table::fread(filename,  header=T, quote="\"", sep=",",na.strings = "-9999",integer64="character")
#   #eur
#   # site<-do.call(rbind,strsplit(basename(filename),'_'))[,2]
#   #ameri
#   site<-do.call(rbind,strsplit(basename(filename),'_'))[,2]
#   kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
#   kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
#   #############################################################################################################
#   ###1. organize the df
#   #############################################################################################################
#   #choose NEE column if null, choose FC (carbon flux column)
#   if(length(grep("FC", names(fluxnet_data), value = TRUE))>=length(grep("NEE", names(fluxnet_data), value = TRUE))){
#     namec<-grep("FC", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$NEE_f<-fluxnet_data[,..namec]
#     fluxnet_data$NEE_fqc=rep(0,length(fluxnet_data[,1])) 
#   }else{
#     namec<-grep("NEE", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$NEE_f<-fluxnet_data[,..namec]
#     fluxnet_data$NEE_fqc=rep(0,length(fluxnet_data[,1])) 
#   }
#   #some NEE_PI report NAs..??
#   if(sum(is.na(fluxnet_data$NEE_f))==length(fluxnet_data$NEE_f)){
#     namec<-grep("FC", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$NEE_f<-fluxnet_data[,..namec]
#   }
#   ##### correct the names
#   ##sw_in
#   if(length(grep("^SW_IN$", names(fluxnet_data), value = TRUE))==0){
#     namesw<-grep("^SW_IN", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$SW_IN<-fluxnet_data[,..namesw]
#     
#   }
#   ##netrad
#   if(length(grep("^NETRAD$", names(fluxnet_data), value = TRUE))==0){
#     namenetr<-grep("^NETRAD", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$NETRAD<-fluxnet_data[,..namenetr]
#     
#   }
#   ##USTAR
#   if(length(grep("^USTAR$", names(fluxnet_data), value = TRUE))==0){
#     nameneustr<-grep("^USTAR", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$USTAR<-fluxnet_data[,..nameneustr]
#     
#   }
#   # ##VPD
#   # if(length(grep("^VPD$", names(fluxnet_data), value = TRUE))==0){
#   # 	namenevpd<-grep("^VPD", names(fluxnet_data), value = TRUE)[1]
#   # 	fluxnet_data$VPD_PI<-fluxnet_data[,..namenevpd]
#   # 	
#   # }
#   #some stattions report ony PPFD and not SW_IN
#   if((length(grep("SW_IN", names(fluxnet_data), value = TRUE))==0 & length(grep("PPFD_IN", names(fluxnet_data), value = TRUE))>=1)){
#     nameppfd<-grep("PPFD_IN", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$SW_IN<-fluxnet_data[,..nameppfd]/(kfFEC*(1 - kalb_vis))
#     fluxnet_data$SW_IN[fluxnet_data$SW_IN<0]<-0
#   }else if(length(grep("SW_IN", names(fluxnet_data), value = TRUE))>=1 & sum(is.na(fluxnet_data$SW_IN))==length(fluxnet_data$SW_IN)&length(grep("PPFD_IN", names(fluxnet_data), value = TRUE))>=1){
#     nameppfd<-grep("PPFD_IN", names(fluxnet_data), value = TRUE)[1]
#     fluxnet_data$SW_IN<-fluxnet_data[,..nameppfd]/(kfFEC*(1 - kalb_vis))
#     fluxnet_data$SW_IN[fluxnet_data$SW_IN<0]<-0
#   }
#   ####get TA average
#   if(length(grep("^TA$", names(fluxnet_data), value = TRUE))==0){
#     nameta<-grep("^TA_", names(fluxnet_data), value = TRUE)
#     fluxnet_data$TA<-rowMeans(fluxnet_data[,..nameta],na.rm=T)
#     
#     
#   }
#   ####get RH average
#   if(length(grep("^RH$", names(fluxnet_data), value = TRUE))==0){
#     namerh<-grep("^RH_", names(fluxnet_data), value = TRUE)
#     fluxnet_data$RH<-rowMeans(fluxnet_data[,..namerh],na.rm=T)
#     fluxnet_data$RH[fluxnet_data$RH>100]<-99.99
#     
#   }
#   #calc vpd hPa if not provided
#   if(length(grep("^VPD", names(fluxnet_data), value = TRUE))==0){
#     fluxnet_data$VPD_PI <- fCalcVPDfromRHandTair(fluxnet_data$RH, fluxnet_data$TA)
#     
#   }else if(length(grep("^VPD", names(fluxnet_data), value = TRUE))>=1){
#     namevpd<-grep("^VPD_", names(fluxnet_data), value = TRUE)
#     fluxnet_data$VPD_PI<-rowMeans(fluxnet_data[,..namevpd],na.rm=T)
#     #fluxnet_data$VPD_PI<-fluxnet_data$VPD
#   }
#   fluxnet_data$VPD_PI[fluxnet_data$VPD_PI<0]<-0.0001
#   
#   #FC, SW_IN, RH, TA, USTAR, L and E
#   
#   
#   #############################################################################################################
#   ###2. build the object required by REddyProc
#   #############################################################################################################	
#   
#   
#   
#   ind<-as.POSIXct(strptime(fluxnet_data$TIMESTAMP_START,format="%Y%m%d%H%M",tz="GMT"))
#   yrs<-unique(format(ind,'%Y'))
#   fluxnet_data$DateTime<-ind
#   if(all(!is.na(start_q),!is.na(start_q))){
#     fluxnet_data <- fluxnet_data %>% filter(TIMESTAMP_START>start_q, TIMESTAMP_END<=end_q)
#   }
#   # fluxnet_data[which(fluxnet_data$NEE_f > 100),'NEE_f'] <- NA
#   
#   EddyProc.C <- sEddyProc$new(site,fluxnet_data, c('NEE_f','NEE_fqc','NETRAD','TA','VPD_PI', 'USTAR','SW_IN'))
#   ## gap fill NEE
#   EddyProc.C$sMDSGapFill('NEE_f', FillAll.b=TRUE)
#   ## get time zone to calc potential rad thus NT
#   if(lon<0){
#     offset = -1 * floor(lon * 24 / 360)
#   }else{
#     offset =  floor(lon * 24 / 360)
#   }
#   
#   EddyProc.C$sSetLocationInfo(LatDeg=lat, LongDeg=lon, TimeZoneHour=offset)
#   EddyProc.C$sEstimateUstarScenarios( nSample = 100L, probs = c(0.05,0.5,0.95))#,seasonFactor = seasonFactor)
#   EddyProc.C$useSeaonsalUStarThresholds()
#   EddyProc.C$sGetUstarScenarios()
#   EddyProc.C$sMDSGapFillUStarScens('NEE_f', FillAll = TRUE)
#   EddyProc.Cc$sFillVPDFromDew()
#   #get GPP NT
#   EddyProc.C$sMRFluxPartition(TempVar = "TA", QFTempVar = "NEE_fqc", RadVar = "SW_IN", parsE0Regression=list(TempRange=1L))
#   EddyProc.C$sMRFluxPartitionUStarScens()
#   #if GPP NT did not work, get GPP_DT
#   
#   if(is.null(EddyProc.C@.xData$sTEMP$GPP_f)){
#     EddyProc.C@.xData$sDATA$Rg_fqc<-EddyProc.C@.xData$sDATA$NEE_fqc
#     EddyProc.C@.xData$sDATA$NEE_fsd<-EddyProc.C@.xData$sTEMP$NEE_f_fsd
#     EddyProc.C$sGLFluxPartition(TempVar = "TA", QFTempVar = "NEE_fqc", RadVar = "SW_IN", VPDVar ='VPD_PI',controlGLPart = partGLControl(nBootUncertainty=0L, isAssociateParmsToMeanOfValids=FALSE , isLasslopPriorsApplied=TRUE, isBoundLowerNEEUncertainty=FALSE , smoothTempSensEstimateAcrossTime=FALSE,useNightimeBasalRespiration=T))
#     EddyProc.C$sGLFluxPartitionUStarScens()
#     fluxnet_data$GPP<-EddyProc.C@.xData$sTEMP$GPP_DT
#     fluxnet_data$RECO<-EddyProc.C@.xData$sTEMP$Reco_DT
#     fluxnet_data$GPP_05 <- EddyProc.C@.xData$sTEMP$GPP_DT_U05
#     fluxnet_data$GPP_95 <- EddyProc.C@.xData$sTEMP$GPP_DT_U95
#     # fluxnet_data$SW_IN_POT<-EddyProc.C@.xData$sTEMP$PotRad
#     fluxnet_data$DateTime <- NULL
#     #write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_DT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
#     write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_DT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
#     
#   }else{
#     fluxnet_data$GPP<-EddyProc.C@.xData$sTEMP$GPP_f
#     fluxnet_data$RECO<-EddyProc.C@.xData$sTEMP$Reco
#     fluxnet_data$GPP_05 <- EddyProc.C@.xData$sTEMP$GPP_U05
#     fluxnet_data$GPP_95 <- EddyProc.C@.xData$sTEMP$GPP_U95
#     # fluxnet_data$SW_IN_POT<-EddyProc.C@.xData$sTEMP$PotRad_NEW
#     fluxnet_data$DateTime <- NULL
#     #write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_NT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
#     write.csv(fluxnet_data,paste0(outpath,'/','AMF_L4_NT_', site, '_HH_', yrs[1],'-',yrs[length(yrs)],'.csv'), row.names=FALSE)
#   }
#   
#   
# }
# 

# 
# 
# 
# checkfc<-function(fluxnet_data){
# 	return(grep("GPP", names(fluxnet_data), value = TRUE))
# }
# 
# 
# 
# checkfc(amerifluxdata[[1]])
# 
# neeavail<-sapply(amerifluxdata,checkfc)
# 
# 
# testf<-grep("FC", names(fluxnet_data), value = TRUE)
# 
# 
# amerifluxdata<-mapply(data.table::fread,filenames.ameriflux[,1],MoreArgs = list( header=T, quote="\"", sep=",",na.strings = "-9999",integer64="character"),SIMPLIFY = F)

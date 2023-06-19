###############################################################################
# get spp at amfx
###############################################################################
get_spp_AMF<-function(site,BADM){
  site_md<-subset(BADM,BADM$SITE_ID==site)
  ###overstory species
  spp_site_o<-subset(site_md,site_md$VARIABLE=='SPP_O')
  spp_cob_site_o<-subset(site_md,site_md$VARIABLE=='SPP_O_PERC')
  spp_date_o<-subset(site_md,site_md$VARIABLE_GROUP=='GRP_SPP_O' & site_md$VARIABLE=='SPP_DATE')
  spp_result_o<-merge(spp_site_o,spp_cob_site_o,by= 'GROUP_ID')
  ### NAs if no date
  if(length(spp_date_o[,1])==0){
    spp_result_o$DATAVALUE<-rep(NA,length(spp_result_o[,1]))
  }else{
    spp_result_o<-merge(spp_result_o,spp_date_o,by= 'GROUP_ID')
  }
  
  ## understory species
  
  spp_site_u<-subset(site_md,site_md$VARIABLE=='SPP_U')
  spp_cob_site_u<-subset(site_md,site_md$VARIABLE=='SPP_U_PERC')
  spp_date_u<-subset(site_md,site_md$VARIABLE_GROUP=='GRP_SPP_U' & site_md$VARIABLE=='SPP_DATE')
  spp_result_u<-merge(spp_site_u,spp_cob_site_u,by= 'GROUP_ID')
  
  ### NAs if no date
  if(length(spp_date_u[,1])==0){
    spp_date_u$DATAVALUE<-rep(NA,length(spp_result_u[,1]))
  }else{
    spp_result_u<-merge(spp_result_u,spp_date_u,by= 'GROUP_ID')
  }
  ##merge all
  spp_result<-data.frame(spp=c(spp_result_o$DATAVALUE.x,spp_result_u$DATAVALUE.x),cover= c(spp_result_o$DATAVALUE.y,spp_result_u$DATAVALUE.y ),
                         date= c(spp_result_o$DATAVALUE,spp_result_u$DATAVALUE ))
  
  spp_result$cover<-as.numeric(as.character(spp_result$cover))
  spp_result
}


###############################################################################
##correct spp names from NRCS
###############################################################################
fix_NRCS_spp<-function(df,NRCS_list){
  sppNRCS_spp<- do.call(rbind,strsplit(as.character(df$spp),' (NRCS plant code)',fixed = T))[,1]
  spp_real<-NRCS_list$Scientific.Name.with.Author[match(sppNRCS_spp,NRCS_plant_code$Symbol)]
  spp_real_alt<-NRCS_list$Scientific.Name.with.Author[match(sppNRCS_spp,NRCS_list$Synonym.Symbol)]
  spp_real[is.na(spp_real)]<-spp_real_alt[is.na(spp_real)]
  ### keep the spp not found
  spp_real[is.na(spp_real)]<-sppNRCS_spp[is.na(spp_real)]
  df$spp<-spp_real
  df <- df[order(df$date,df$cover,decreasing = T),]
  df
}
###############################################################################
#get spps NEON sites
###############################################################################
get_text_pars<-function(x){
  ###neon sites
  neon_sites<-regmatches(x, gregexpr("NEON*.*", x, perl=T))
  neon_sites<-lapply(neon_sites,function(x){if(length(x)==0){NA}else{x}})
  neon_sites<-do.call(c,neon_sites)
  result<-regmatches(neon_sites, gregexpr("(?<=\\().*?(?=\\))", neon_sites, perl=T))
  result<-lapply(result,function(x){if(length(x)==0){NA}else{x}})
  result<-do.call(c,result)
  # if(length(result)==0){
  # 	result=NA
  # }
  return(result)
}
###############################################################################
#function ge spp cover NEON
###############################################################################
get_spp_NEON<-function(filename){
  spp_df<-read.table(file = filename, header = TRUE, sep = ",", fileEncoding = "windows-1252", quote = "\"", stringsAsFactors = FALSE, comment.char = "", na.strings = "")
  #### add non vegetation cover
  spp_df$scientificName[is.na(spp_df$scientificName)]<-'non_photo(bare soil_dead wood,litter)'
  spp_df$endDate<-as.Date(spp_df$endDate)
  spp_cov<-aggregate(spp_df$percentCover, by = list(spp_df$scientificName), mean,na.rm=T)
  dates<-aggregate(spp_df$endDate, by = list(spp_df$scientificName), median,na.rm=T)
  spp_cov<-merge(spp_cov,dates, by='Group.1')
  names(spp_cov)<-c('spp','cover','date')
  spp_cov$cover<-(spp_cov$cover/sum(spp_cov$cover,na.rm=T))*100
  spp_cov <- spp_cov[order(spp_cov$cover,decreasing = T),]
  
  spp_cov
  
}
###############################################################################
##fix species name
###############################################################################

library(Taxonstand)
fix_name<-function(df){
  
  fixed_names <- TPL(df$spp, corr = TRUE)
  df$spp<-fixed_names$Taxon
  df
  
}
###############################################################################
#### read try db
###############################################################################

read_try<-function(spp_name,trait, dbase){
  names_db<-names(dbase)
  ##get most similar name
  spp_name<-unlist(strsplit(spp_name, " "))[1:2]
  spp_name<-paste(spp_name[1],spp_name[2])
  # Find similar (not exact) matches in the target strings using the Jaro-Winkler distance
  distances <- stringdist(spp_name, dbase$SpeciesName, method = "jw")
  # Set a threshold for the minimum distance to consider a match
  threshold <- 0.25
  ###subset by spp
  dbase<-subset(dbase,distances < threshold)
  # spp_name_sim<-(agrepl(spp_name,dbase$SpeciesName,ignore.case = TRUE, max.distance = 0.25))
  # dbase<-subset(dbase ,spp_name_sim)
  #dbase<-subset(dbase,dbase$SpeciesName %in% spp_name_sim)
  #availtra<-unique(dbase$AccSpeciesName)
  ##subset by trait
  dbase<-subset(dbase ,dbase$TraitName==trait)
  dbase$OrigValueStr<-as.numeric(dbase$OrigValueStr)
  dbase<-subset(dbase ,dbase$OrigValueStr==max(dbase$OrigValueStr,na.rm = T))
  if (length(dbase[,1])>1){
    dbase<-dbase[1,]
  }else if (length(dbase[,1])==0){
    dbase[1,]<-rep(NA,28)
  }
  
  dbase
}
###############################################################################
#### match photo traits
###############################################################################
library(stringdist)
find_QY_TRY<-function(df,dbase = all_photo_traits){
  
  
  #get QY
  QY<-mapply(read_try,df$spp,MoreArgs = list(trait='Leaf photosynthesis light use efficiency (LUE)',dbase = dbase),SIMPLIFY = F)
  QY<-do.call(rbind,QY)
  #get QY2
  QY2<-mapply(read_try,df$spp,MoreArgs = list(trait='Leaf photosynthesis quantum yield (QY; corresponding to photosynthetic efficiency)',dbase = dbase),SIMPLIFY = F)
  QY2<-do.call(rbind,QY2)
  #get QY3
  QY3<-mapply(read_try,df$spp,MoreArgs = list(trait='Leaf photosynthesis: quantum efficiency of photosystem II (maximum potential; FvFm)',dbase = dbase),SIMPLIFY = F)
  QY3<-do.call(rbind,QY3)
  
  df$phi_0.spp<-as.numeric(QY$OrigValueStr)
  df$PHI_II.spp<-as.numeric(QY2$OrigValueStr)
  df$PHI_FvFm.spp<-as.numeric(QY3$OrigValueStr)
  df
  
}
# evaluation sim plot

plot.eval.dens<-function(sim,obs, factor=NULL,cex.stats=1.0, ...){
	# testing
	# sim=wn_sim_well;obs=wn_obs_well$sm;igbp=wn_obs_well$IGBP
	# end testing
	sim<-as.numeric(sim)
	obs<-as.numeric(obs)
	#get stats
	stats<-hydroGOF::gof(as.numeric(sim),as.numeric(obs),digits=3, na.rm = TRUE)
	#stats<-c(stats[17,],stats[4,],stats[9,])
	if(!is.null(factor)){
		factor=as.factor(factor)
		sim_group<-aggregate(sim,by=list(factor),FUN=mean,na.rm=TRUE)
		sim_group_sd<-aggregate(sim,by=list(factor),FUN=sd,na.rm=TRUE)
		obs_group<-aggregate(obs,by=list(factor),FUN=mean,na.rm=TRUE)
		obs_group_sd<-aggregate(obs,by=list(factor),FUN=sd,na.rm=TRUE)
		types<-obs_group$Group.1
	}

	#get lm
	fit<-lm(as.numeric(obs)~as.numeric(sim))
	subtitle <- bquote( italic(R)^2 == .(round(stats[17], digits = 2)) ~~
		RMSE == .(round(stats[4], digits = 2)) ~~
		bias == .(round(stats[1], digits = 2)) ~~
		slope == .(round(fit$coefficients[2], digits = 2)) ~~
		italic(N) == .(length(sim)) )	
	
	symbs<-c(15:18,21:25,3:4,9:10,12)
	##plot
	LSD::heatscatter(sim,obs,colpal=c('grey80', 'blue', 'green', 'yellow', 'red'),main='',bty='l',...)
	##add regression lines to the plot
	abline(0,1,lwd=1)
	abline(coef(fit)[1:2],lty=2,lwd=1,col="red")
	##add the stats to the plot
	mtext(subtitle,side=3, line=1, adj=0, cex=cex.stats,outer=F)
	if(!is.null(factor)){
		arrows(sim_group$x-sim_group_sd$x,obs_group$x,sim_group$x+sim_group_sd$x,obs_group$x , length=0.05, angle=90, code=3,col=1)
		arrows(sim_group$x,obs_group$x-obs_group_sd$x,sim_group$x,obs_group$x+obs_group_sd$x , length=0.05, angle=90, code=3,col=1)
		points(sim_group$x,obs_group$x,pch=symbs[as.integer(types)],cex=1.7,col='black',bg='grey')
		legend("bottomright",pch=symbs[as.integer(types)],col=rep('black',length(types)),legend=types,bty="n",pt.bg=rep('grey',length(types)),cex=1.3,pt.cex=1.5)
	}
	
}

	#harvestFun <- function(t, y, parms){
	#	B<- y
	#	C<- harvest[harvest$time==t,]$value
	#	x<- unlist(ifelse(abs(C)<B, (B+C), 0.000000000000000000000001))
	#	B <- as.numeric(x)[1]
	#	return(c(B))}	
		
	harvestFun <- function(t, y, parms){
		B<- y
		F<- as.numeric(parms['F'])
		B <-as.numeric(B - ifelse(round(t- trunc(t),1) == 0.2, B*F*0.86, B*F*0.14))
		return(c(B))}	

	semidiscrete <- function(time, y, parms){
		B<- ifelse(y<0, 0.0001, y)
		r<- parms['r']
		Bmax<- parms['Bmax']	
		p<- parms['p']
		m<- parms['m']
		if(m==1){dB<-r*B}
		if(m==2){dB<-r*B*((Bmax-B)/Bmax)}
		if(m==3){dB<-r*B*(log(Bmax/B))}
		if(m==4){dB<-(r/p)*B*(1-(B/Bmax)^p)}
		if(m==5){dB<-r*B*((Bmax-B)/Bmax)}
		return(list(c(dB)))
		}	
	#continuous <- function(time, y, parms){ 
	#	B<- ifelse(y<0, 0.0001, y)
	#	r<- parms['r']
	#	Bmax<- parms['Bmax']	
	#	p<- parms['p']
	#	m<- parms['m']
	#	C <- -1*harvestConFunction(time)
	#	if(m==1){dB<- r*B - C}
	#	if(m==2){dB<- r*B*((Bmax-B)/Bmax) - C}
	#	if(m==3){dB<- r*B*(log(Bmax/B)) - C}
	##	if(m==4){dB<- (r/p)*B*(1-(B/Bmax)^p) - C}
		#if(m==5){dB<- r*B*((Bmax-B)/Bmax) - C}
	#	dB<- ifelse(C<B,dB, -1*B)
	#	return(list(c(dB),Y=C))
	#	}	

	continuous <- function(time, y, parms){ 
		B<- y
		r<- parms['r']
		F<- parms['F']
		Bmax<- parms['Bmax']	
		p<- parms['p']
		m<- parms['m']
		C <- F*B
		if(m==1){dB<- r*B - C}
		if(m==2){dB<- r*B*((Bmax-B)/Bmax) - C}
		if(m==3){dB<- r*B*(log(Bmax/B)) - C}
		if(m==4){dB<- (r/p)*B*(1-(B/Bmax)^p) - C}
		if(m==5){dB<- r*B*((Bmax-B)/Bmax) - C}
		dB<- ifelse(C<B,dB, -1*B)
		return(list(c(dB),Y=C))
		}	

	#----------------------------------------------------------------------------------------------
	# LIKLIHOOD FUNCTION FOR CONTIONUOUS MODEL
	#----------------------------------------------------------------------------------------------	
	lnLikFunc_con<- function(initParms,data){
		r<- initParms['r']
		p<- initParms['p']
		Bmax<- initParms['Bmax']
		F<- parms['F']
		parms<-c(Bmax, p, r, m = dat$model,F)		
		initialValues<- c(B= as.numeric(initParms['Bmax']))
		sigma<- initParms['sigma']
		obs<- dat$obs
		modelout  <-ode(
			y 		= initialValues, 
			times 	= times, 
			func 	= continuous, 
			parms 	= parms,
			events 	= NULL,
			method 	= "lsoda")
		modelout<- data.frame(modelout)
		predictedbiomass<- approxfun(modelout$time, modelout$B, method="linear", rule=2)	
		pred<-predictedbiomass(sampleTimes)
		pred[pred<=0]<- 0.000000000000000000000001
		pred[pred=="NaN"]<- 0.000000000000000000000001
		resids<- log(pred/obs)
		lnLik= sum(dnorm(resids,0,sigma,log=T))
		return(-1*lnLik)}	
	#----------------------------------------------------------------------------------------------
	# LIKLIHOOD FUNCTION FOR SEMI-DISCRETE MODEL
	#----------------------------------------------------------------------------------------------		
	lnLikFunc_semiD<- function(initParms, data){
		r<- initParms['r']
		p<- initParms['p']
		Bmax<- initParms['Bmax']
		F<- parms['F']
		parms<-c(Bmax, p, r, m = dat$model,F)	
		initialValues<- c(B= as.numeric(initParms['Bmax']))
		sigma<- initParms['sigma']
		obs<- dat$obs
		modelout  <-ode(
				y 		= initialValues, 
				times 	= times, 
				func 	= semidiscrete, 
				parms 	= parms,
				events 	= list(func= harvestFun, time=evtimes),
				method 	= "lsoda")
		modelout<- data.frame(modelout)
		predictedbiomass<- approxfun(modelout$time, modelout$B, method="linear", rule=2)	
		pred<-predictedbiomass(sampleTimes)
		pred[pred<=0]<- 0.000000000000000000000001
		pred[pred=="NaN"]<- 0.000000000000000000000001
		resids<- log(pred/obs)
		lnLik= sum(dnorm(resids,0,sigma,log=T))
		return(-1*lnLik)}	
	#----------------------------------------------------------------------------------------------
	# LIKLIHOOD FUNCTION FOR CONTIONUOUS ASPIC MODEL
	#----------------------------------------------------------------------------------------------	
	lnLikFunc_con_aspic<- function(initParms,data){
		r<- initParms['r']
		#p<- initParms['p']
		Bmax<- initParms['Bmax']
		initialValues<- c(B=as.numeric(initParms['Bmax']))
		q<- initParms['q']
		F<- parms['F']
		parms<-c(Bmax, r,m = dat$model,F)	
		sigma<- initParms['sigma']
		sigma_cpue<- initParms['sigma_cpue']
		obs		<- dat$obs
		obs_cpue<-dat$cpue
		modelout  <-ode(
			y 		= initialValues, 
			times 	= times, 
			func 	= continuous, 
			parms 	= parms,
			events 	= NULL,
			method 	= "lsoda")
		modelout<- data.frame(modelout)
		predictedbiomass<- approxfun(modelout$time, modelout$B, method="linear", rule=2)	
		pred<-predictedbiomass(sampleTimes)
		cpue_pred<- q*predictedbiomass(sampleTimes_cpue)
		pred[pred<=0]<- 0.000000000000000000000001
		pred[pred=="NaN"]<- 0.000000000000000000000001
		resids<- log(pred/obs)
		cpue_pred[cpue_pred<=0]<- 0.000000000000000000000001
		cpue_pred[cpue_pred=="NaN"]<- 0.000000000000000000000001
		resids_cpue<- log(cpue_pred/obs_cpue)
		lnLik= sum(dnorm(resids,0,sigma,log=T))+ sum(dnorm(resids_cpue,0,sigma_cpue,log=T))
		return(-1*lnLik)}	
	#----------------------------------------------------------------------------------------------
	# LIKLIHOOD FUNCTION FOR SEMI-DISCRETE ASPIC MODEL
	#----------------------------------------------------------------------------------------------		
	lnLikFunc_semiD_aspic<- function(initParms, data){
		r<- initParms['r']
		#p<- initParms['p']
		Bmax<- initParms['Bmax']
		initialValues<- c(B=as.numeric(initParms['Bmax']))
		q<- initParms['q']
		F<- parms['F']
		parms<-c(Bmax, r,m = dat$model,F)
		sigma<- initParms['sigma']
		sigma_cpue<- initParms['sigma_cpue']
		obs		<- dat$obs
		obs_cpue<- dat$cpue
		modelout  <-ode(
				y 		= initialValues, 
				times 	= times, 
				func 	= semidiscrete, 
				parms 	= parms,
				events 	= list(func= harvestFun, time=evtimes),
				method 	= "lsoda")
		modelout<- data.frame(modelout)
		predictedbiomass<- approxfun(modelout$time, modelout$B, method="linear", rule=2)	
		pred<-predictedbiomass(sampleTimes)
		cpue_pred<- q*predictedbiomass(sampleTimes_cpue)
		pred[pred<=0|pred=="NaN"]<- 0.000000000000000000000001
		#pred[pred=="NaN"]<- 0.000000000000000000000001
		resids<- log(pred/obs)
		cpue_pred[cpue_pred<=0|cpue_pred=="NaN"]<- 0.000000000000000000000001
		#cpue_pred[cpue_pred=="NaN"]<- 0.000000000000000000000001
		resids_cpue<- log(cpue_pred/obs_cpue)
		lnLik= sum(dnorm(resids,0,sigma,log=T))+ sum(dnorm(resids_cpue,0,sigma_cpue,log=T))
		return(-1*lnLik)}	
	
	
	
	

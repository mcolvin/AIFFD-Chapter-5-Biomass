#===========================================================================================================================
# THIS IS A FUNCTION TO CALCULATE THE MASS BALANCE
# SOLUTION OF AN ECOPATH MODEL
# 
#
# 27 NOVEMBER 2010 COLVIN.MIKE@GMAIL.COM
#===========================================================================================================================
	basicInputs<- dat
	rownames(basicInputs)<- basicInputs$functionalGroup
	nFG<- nrow(dc)
	nFG <- nFG-1 # get detritus out of there
	nConsumers<-nrow(basicInputs[basicInputs$type==2,])
	matrixNames<- basicInputs$functionalGroup[-nrow(basicInputs)]
	#==============================================================================================
	# A ~= b: APPROXIMATE EQUALITIES
	#==============================================================================================
		p	<- (basicInputs$biomass*basicInputs$pb)[-nrow(basicInputs)]; 
		q	<- (basicInputs$biomass*basicInputs$qb)[1:nConsumers] # vector of consumptions for consumers	
		ba	<- (basicInputs$ba)[-nrow(basicInputs)]
		y	<- (basicInputs$comHarvest + basicInputs$recHarvest)[-nrow(basicInputs)]
		ex	<-( basicInputs$export)[-nrow(basicInputs)]
		si	<- (basicInputs$stocking + basicInputs$import)[-nrow(basicInputs)]
		#Q<- as.matrix(DC);Q[Q<0]<-1
		#Y[Y<0]<-1;Ex[Ex<0]<-1
		#A<- adiag(P,Q, BA, Y,Ex,SI)

		names<- c(paste(matrixNames,"ee",sep="-"),
			paste(matrixNames[1:nConsumers],"consumption",sep="-"),
			paste(matrixNames,"ba",sep="-"),
			paste(matrixNames,"y",sep="-"),
			paste(matrixNames,"export",sep="-"),
			paste(matrixNames,"import",sep="-"))
		b	<- c(p,q,ba,y,ex,si)			
		A<- diag(1, nrow=length(b), ncol=length(b))
	# RETAIN NON ZERO VALUES
		A<- A[b != 0,]
		b	<- b[b != 0]
		
	#==============================================================================================	
	# E = f: MATRICES FOR EQUALITIES 
	#==============================================================================================
		P	<- diag( 1,nrow=nFG, ncol=nFG) # a base matrix of diagonal 1s 		
		DC	<- as.matrix(dc[c(1:nFG),c(1:nConsumers)]*-1)
		BA	<- diag( 1,nrow=nFG, ncol=nFG) 
		Y	<- diag(-1,nrow=nFG, ncol=nFG)
		Ex	<- diag(-1,nrow=nFG, ncol=nFG)
		SI	<- diag( 1,nrow=nFG, ncol=nFG) # stocking and import from ventura marsh
	# E MATRIX
		E	<- cbind(P, DC, BA, Y, Ex, SI)
	# FILL ADDITIONAL EQUALITIES (E.G., Y = 0)
		XXX<- diag(1, nrow=length(b), ncol=ncol(E))
		XXX<- XXX[b==0,]
		E <- rbind(E, XXX)
		colnames(E)<- names
	# EQUALITY RATES
		f<- rep(0, nrow(E))

	#==============================================================================================
	# inequalities (constrain some to be positive)
	# constrain ecotrophic efficiencies to be between 0 and 1
	#==============================================================================================
	# ECOTROPHIC EFFICIENCIES
		eelower<- diag(1,nrow=nFG, ncol= nFG)
		eeupper<- diag(-1,nrow=nFG, ncol= nFG)
		EE<- rbind(eelower, eeupper)
		ee<- c(rep(0,nFG),rep(-1,nFG)) 
	# constrain consumption to be greater than 0
		QLOW<- diag(1, nrow=nConsumers, ncol= nConsumers)
		qlow<- rep(0.01, nConsumers)
	# BA AS A PLACE HOLDER
		BAHOLD<- matrix(0,nrow=nFG, ncol=nFG) 
		bahold<- rep(0,nFG)
	# Y
		YCONSTRAINT<- diag(1,nrow=nFG, ncol=nFG) 
		yconstraint<- rep(0,nFG)	
	# IMPORTS
		IMPORTS<- diag(1,nrow=nFG, ncol=nFG) 
		imports<- rep(0,nFG)
	# EXPORTS
		EXPORTS<- diag(1,nrow=nFG, ncol=nFG) 
		exports<- rep(0,nFG)
	
	G <- adiag(EE,QLOW, BAHOLD, YCONSTRAINT, EXPORTS, IMPORTS)
	h <- c(ee, qlow, bahold, yconstraint, exports, imports)

	#==============================================================================================
	# SOLVE THE LINEAR INVERSE MODEL	
	#==============================================================================================
	res<- lsei( E=E, F=f, G=G, H=h)# A=A, B=b,
	origParms	<- c(p,q,ba,y,ex,si)
	out<- cbind(res$X, origParms)
	plot(res$X~origParms)
	ee<- res$X[c(1:nFG)]
	q<- res$X[c(nFG+1:nFG+29)]
	matrix(res$X, nrow=34, byrow=FALSE)
	
	xs<- xsample(A=A, B=b, E=E, F=f, G=G, H=h,sdB=NULL)
	pairs(xs$X)
	xranges(E=E, F=f, G=G, H=h, central=TRUE)

	error <- qt(0.975,df=1)*s/sqrt(n)
	s<- error/qt(0.975,df=1) # get sd from error s=sd

	
	
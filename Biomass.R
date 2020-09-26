

box1<- data.frame(
    year=c(1985,1986,1987,1988,1989,1990,1991,1992, 
        1993,1994,1995,1996,1997,1998,1999,2000,2001)
    effort=c(825,1008,1411,1828,2351,2074,1877,1566,1139, 
        893, 1029,727,658,953,1012,1203,1034)
    catch=c(90000,113300,155860,181128,198584,198395,139040, 
        109969,71896,59314,62300,65343,76990,88606,118016,108250,108674)

box1$cpue<- box1$catch/box1$effort


B0 <- 800000                                  # initial value of B0
K <- 1000000                                  # initial value of K
q <- 0.0001                                   # initial value of q
r <- 0.17                                     # initial value of r
pars <- c(B0,K,q,r)                           # put all parameters into one vector
names(pars) <- c("B0","K","q","r")            # name the parameters

SPsse <- function(par,B,CPE,SSE.only=TRUE)  { ## This repeats Excel calculations
  n <- length(B)                              # get number of years of data
  B0<- par["B0"]                             # isolate B0 parameter
  K <- par["K"]                               # isolate K parameter
  q <- par["q"]                               # isolate q parameter
  r <- par["r"]                               # isolate r parameter
  predB <- numeric(n)
  predB[1] <- B0
  for (i in 2:n) predB[i] <- predB[i-1]+r*predB[i-1]*(1-predB[i-1]/K)-B[i-1]
  predCPE <- q*predB
  sse <- sum((CPE-predCPE)^2)
  if (SSE.only) sse
    else list(sse=sse,predB=predB,predCPE=predCPE)
}
( res1 <- SPsse(pars,d7$catch,d7$cpe) )

pars.scale <- c(1e5,1e6,1e-4,1e-1)            # set rescale values for parameters
SPoptim1 <- optim(pars,SPsse,
    control=list(maxit=10000,parscale=pars.scale),
    B=d7$catch,
    CPE=d7$cpe)
SPoptim1$value                                # "minimum" SSE
SPoptim1$counts[1]                            # number of iterations
SPoptim1$convergence                          # a "0" means completed successfully


presPlot()
plot(catch~year,d7,ylab="Catch (kg)",xlab="Year",type='b')

res3 <- SPsse(SPoptim1$par,d7$catch,d7$cpe,SSE.only=FALSE)
str(res3)
plot(cpe~year,data=d7,pch=19,xlab="Year",ylab="Catch per unit effort")
lines(d7$year,res3$predCPE,lwd=2,col="red")


plot(res3$predB~d7$year,pch=19,xlab="Year",ylab="Predicted biomass (kg)",type='b')


r=0.41
K=1158451
B0= 706743
F<-seq(0,0.5,by=0.01)
tdy<-list()
out<-list()
for(j in 1:length(F))
    {
    B<-rep(0,100)
    H<-rep(0,100)
    B[1]<- B0
    for(i in 2:100)
        {
        B[i] <- B[i-1]+r*B[i-1]*(1-B[i-1]/K)-F[j]*B[i-1]
        H[i-1]<- F[j]*B[i-1]
        }
    tdy[[j]]<- data.frame(year=c(1:100),F=F[j],B=B,H=H)
    out[[j]]<- data.frame(F=F[j],B=B[99],H=H[99])
    }

tdy=do.call("rbind",tdy)
out=do.call("rbind",out)

lattice::xyplot(B~year,tdy, groups=F,type='l',
    xlab="Year",ylab="Biomass (kg)")

out[which.max(out$H),]
presPlot()
plot(H~F,out,type='l',ylab="Yield (kg)", xlab="Fishing mortality")

## F0.1
slpatorigin<-((out$H[2]-out$H[1])/0.01)
slpatorigin01<-((out$H[2]-out$H[1])/0.01)*0.1
presPlot()
plot(H~F,out,type='l',ylab="Yield (kg)", xlab="Fishing mortality")
abline(0,slpatorigin,lty=2)
abline(0,slpatorigin01,lty=2)

out$F[which.max(slpatorigin01-abs((out$H[-1]-out$H[-length(out$H)])/0.01))]





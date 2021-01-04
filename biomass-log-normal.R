

# Biomass estimates

Does it matter if weight is normally distributed when using mean weight
to estimate biomass as $N\cdot \bar{W}$?


First we need a resaonble population to evaluate. 

The total instantaneous mortality ($Z$) for this population is `r Z<-0.7274;print(Z)`, 
with natural mortaltity ($M$) of `r M<-0.6009;print(M)`


```{r}
library(data.table)
age<-c(1:4)
N<-10000
for(i in 2:length(age))
    {
    N<-c(N,N[i-1]*exp(-1*Z*i))
    }
dat<-data.table(A=age,N=round(N,0))
plot(N~A,dat)
vbgf<- function(linf, t0, k, age){linf * (1 - exp(-k * (age-t0))) }
dat<-dat[rep(seq(1, nrow(dat)), dat$N),]
dat[,spp:="FLMB"]
dat[,L:=vbgf(linf=550,t0=0,k=0.8,A)]#cm
dat$L<- round(dat$L*rnorm(nrow(dat),1,0.05),2)
plot(L~A,dat)

dat[,W:=10^(-5.528+3.272*log10(L))]
dat$W<- round(dat$W*rnorm(nrow(dat),1,0.1),0)
plot(W~L,dat)
dat<-as.data.table(dat)
q<- rep(0.01,length(age))
q<- c(0.01,0.02,0.025,0.01)
dat[,captured:=rbinom(.N,1,0.01)]
catch<- subset(dat,captured==1)
biomass<-sum(catch$W)
biomass_est<-nrow(catch)*mean(catch$W)
hist(catch$W)







```

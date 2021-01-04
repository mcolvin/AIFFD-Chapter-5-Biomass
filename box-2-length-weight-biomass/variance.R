


```{r,eval=FALSE}
x<-data.frame(W=rnorm(50,1.2,0.1),group=rep(1,50))
N<-5000
```

## Fitting the length weight relationship

```{r}
y<- lm(W~group,x)
mean(x$W)
sqrt(var(x$W)/50)# standard error same as se from predict
pp<-predict(y,newdata=data.frame(group=rep(1,N)),se.fit = TRUE)
sum(pp$se.fit)^2# variance is the sum of predicted ses squared
# variance for biomass 
N^2*(var(x$W)/50)
sum(pp)
```


x<-data.frame(W=rnorm(50,1.2,0.1),group=rep(1,50))
N<-5000
y<- lm(W~group,x)
mean(x$W)
sqrt(var(x$W)/50)# standard error same as se from predict
pp<-predict(y,newdata=data.frame(group=rep(1,N)),se.fit = TRUE,interval = "confidence")
sum(pp$se.fit)^2# variance is the sum of predicted ses squared
# variance for biomass 
N^2*(var(x$W)/50)
sum(pp)
```
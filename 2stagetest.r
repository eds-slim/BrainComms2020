set.seed(1)

n <- 30 #number of subjects
t <- c(0,1,3,12) # time points
b <- 1 # coefficient beta
c <- -2 # coefficient gamma

x0 <- rnorm(n,10,1) # random intercepts

x <- unlist(lapply(t, function(t){x0+b*t}))
y <- c*x

df<-data.frame(x=x+rnorm(4*n,0,1), y=y+rnorm(4*n,0,1), t=rep(t,each=n), ID=rep(seq.int(1,n),4))

require(ggplot2)
ggplot(df,aes(x=x,y=y,col=as.factor(t)))+
  geom_point()+
  geom_smooth(method='lm', se = F)


fm <- lme(y ~ x+t, data = df, random = ~ 1 | ID)
summary(fm)

fm <- lme(x ~ t, data = df, random = ~ 1 | ID)
summary(lm(y ~ 1 + fitted(fm), data = df))

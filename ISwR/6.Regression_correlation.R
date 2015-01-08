### 6.1 Simple linear regression
attach(thuesen)
lm(short.velocity~blood.glucose)
summary(lm(short.velocity~blood.glucose)) #extractor function

plot(blood.glucose,short.velocity)
abline(lm(short.velocity~blood.glucose))

### Extraction
lm.velo=lm(short.velocity~blood.glucose)
fitted(lm.velo)
resid(lm.velo)

plot(blood.glucose,short.velocity)
lines(blood.glucose[!is.na(short.velocity)],fitted(lm.velo))

cc=complete.cases(thuesen)
thuesen[cc,]

segments(blood.glucose,fitted(lm.velo),blood.glucose,short.velocity)
qqnorm(resid(lm.velo))


### Prediction and confidence
predict(lm.velo,int="c")
predict(lm.velo,int="p")

pred.frame=data.frame(blood.glucose=4:20)
pp=predict(lm.velo,int="p",newdata=pred.frame)
pc=predict(lm.velo,int="c",newdata=pred.frame)
plot(blood.glucose,short.velocity,ylim=range(short.velocity,pp,na.rm=T))
pred.gluc=pred.frame$blood.glucose
matlines(pred.gluc,pc,lty=c(1,2,2),col="black")
matlines(pred.gluc,pp,lty=c(1,3,3),col="black")


### Correlations
# Pearson
cor(blood.glucose,short.velocity,use="complete.obs")
cor.test(blood.glucose,short.velocity)

#Spearman
cor.test(blood.glucose,short.velocity,method="spearman")

#Kendall
cor.test(blood.glucose,short.velocity,method="kendall")

### Exercises
#6.1
library(ISwR)
attach(rmr)
plot(body.weight,metabolic.rate)
lm.weight=lm(metabolic.rate~body.weight)
abline(lm(metabolic.rate~body.weight))
fit=fitted(lm.weight)
summary(fit)
summary(lm.weight)
811.2267+7.0595*70
predict(lm.weight,newdata=data.frame(body.weight=70))
confint(lm.weight)
detach(rmr)

#6.2
attach(juul)
head(juul)
age.25=juul[age>25,]
lm(sqrt(igf1)~age,data=age.25)

summary(lm(sqrt(igf1)~age,data=juul,subset=age>25))
detach(juul)

#6.3
attach(malaria)
head(malaria)
summary(lm(log(ab)~age),data=malaria)
plot(log(age),ab)
abline(lm(log(ab)~age,data=malaria))

#6.4
x=rnorm(n=5000,mean=0,sd=1)
y=rnorm(5000,0.95*x,sqrt(1-0.95^2))
cor.test(y,x,method="spearman")
cor.test(y,x,method="kendall")
cor(y,x)
plot(x,y)

# destroy of the correlation
x=c(x,1000)
y=c(y,-1000)
plot(x,y)
cor.test(y,x,method="spearman") ##rho=0.944
cor.test(y,x,method="kendall") ##tau=0.796
cor(y,x) ##t=-0.99!!!

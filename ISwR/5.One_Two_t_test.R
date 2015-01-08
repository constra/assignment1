### One sample test
daily.intake = c(5260,5470,5640,6180,6390,6515,6805,7515,7515,8230,8770)
## Emperical statistics
mean(daily.intake)
sd(daily.intake)
quantile(daily.intake)

## t-test
t.test(daily.intake,mu=7725)
# one-side tests
t.test(daily.intake,mu=7725,alternative="greater")
t.test(daily.intake,mu=7725,alternative="less")
# alternative confidence level
t.test(daily.intake,mu=7725,conf.level=0.99)

## Wilcoxon test
wilcox.test(daily.intake,mu=7725)

### Two sample test
library(ISwR)
attach(energy)

## t-test
t.test(expend~stature)
t.test(expend~stature,var.equal=T) ## assume that variance is the same

detach(energy)
t.test(expend~stature,data=energy)

### Comparison of variances
attach(energy)
var.test(expend~stature)

### Wilcoxon test, two samples
wilcox.test(expend~stature)

### The paired t test
attach(intake)
intake
post-pre
t.test(pre,post,paired=T)
#t.test(pre,post)


### The paired wilcoxon t test
wilcox.test(pre,post,paired=T)
#wilcox.test(pre,post)

### Exercises
#5.1
react
qqnorm(react)
t.test(react,mu=0)
#5.2
vitcap
attach(vitcap)
t.test(vital.capacity~group,conf.level=0.99)
boxplot(age~group) # the cause of misleading result

#5.3
wilcox.test(react)
wilcox.test(vital.capacity~group,conf.level=0.99)

#5.4
intake
plot((pre+post)/2,post-pre,ylim=range(0,post-pre))
abline(h=0)
hist(pre-post)

#5.5
shapiro.test(react)

#5.6
ashina
attach(ashina)
t.test(vas.active,vas.plac,paired=T) #No period effect
t.test((vas.active-vas.plac)[grp==1],(vas.plac-vas.active)[grp==2]) #With period effect

#5.7
t.test(rnorm(25))
t.test(rnorm(25))$p.value

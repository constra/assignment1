library(SDSFoundations)
bull<-BullRiders

plot(bull$Top10,bull$RidePer,xlab="Places Top 10",ylab='How often stayed on the bull',main='Top 10 vs. how often')
abline(lm(bull$RidePer~bull$Top10))
cor(bull$Top10,bull$RidePer)
# Visualize and describe the first variable of interest 
hist(bull$RidePer)
fivenum(bull$RidePer)
mean(bull$RidePer)
sd(bull$RidePer)

# Visualize and describe the second variable of interest 
hist(bull$Top10)
fivenum(bull$Top10)
mean(bull$Top10)
sd(bull$Top10)

# Create a scatterplot
plot(bull$RidePer,bull$Top10)

# Add line of best fit
abline(lm(bull$Top10~bull$RidePer))

# Calculate the correlation coefficient
cor(bull$RidePer,bull$Top10)

# Create a correlation matrix  
vars <- c("Top10", "RidePer")
cor(bull[,vars])
which(bull$Top10==5 & bull$RidePer==.53)


#LAB
hist(bull$Earnings)
vars<-c('Earnings','RidePer','CupPoints')
cor(bull[,vars])
plot(bull$RidePer,bull$Earnings,xlab='Ride per 8sec',ylab='Earnings',main='correlation of Earnings vs Ride')
plot(bull$CupPoints,bull$Earnings,xlab='Cup points',ylab='Earnings',main='correlation of Earnings vs Cup points')
median(bull$Earnings)
max(bull$Earnings)
# identify specific case
which(bull$Earnings == max(bull$Earnings))
#Subset the data
nooutlier <-bull[-1,]
cor(nooutlier[,vars])


###Problem set
RidePerEvents<-bull$RidesPerEvent <- bull$Rides/bull$Events
hist(RidePerEvents)
min(RidePerEvents)
median(RidePerEvents)
max(RidePerEvents)
plot(bull$RidesPerEvent,bull$Place)
abline(lm(bull$Place~bull$RidesPerEvent))

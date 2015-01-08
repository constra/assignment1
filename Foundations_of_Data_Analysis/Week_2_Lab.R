library(SDSFoundations)

animaldata<-AnimalData

# Table showing adult cat and dogs
adult.animals<-animaldata[animaldata$Age.Intake>=1,]
table(adult.animals)
adult.dog<-adult.animals[adult.animals$Animal.Type=="Dog",]
adult.cat<-adult.animals[adult.animals$Animal.Type=="Cat",]

#distribution of the weight

hist(adult.dog$Weight)

hist(adult.cat$Weight)

mean(adult.cat$Weight)
sd(adult.cat$Weight)

#Z-score of the 13 pound aduilt cat
cat.weight.mean<-mean(adult.cat$Weight)
cat.weight.sd<-sd(adult.cat$Weight)
cat.weight.z.13<-(13-cat.weight.mean)/cat.weight.sd

1-pnorm(cat.weight.z.13)

#Dog weight
fivenum(adult.dog$Weight)

#Dog statistics
Dog.data<-animaldata[animaldata$Animal.Type=="Dog",]
table(Dog.data$Intake.Type)

sum(Dog.data$Intake.Type=="Owner Surrender")/length(Dog.data$Intake.Type)
dog.surrender<-Dog.data[Dog.data$Intake.Type=="Owner Surrender",]
dim(dog.surrender[dog.surrender$Outcome.Type=="Return to Owner"])
dog.return<-dog.surrender[dog.surrender$Outcome.Type=="Return to Owner",]
dog.return$Days.Shelter

#Problem set
Dogs<-animaldata[which(animaldata$Animal.Type=="Dog"),]
table(Dogs$Intake.Type)
prop<-sum(Dogs$Intake.Type=="Owner Surrender")/nrow(Dogs)
prop
owner.surrender<-Dogs[which(Dogs$Intake.Type=="Owner Surrender"),]
table(owner.surrender$Outcome.Type)
returned<-owner.surrender[which(owner.surrender$Outcome.Type=="Return to Owner"),]
mean(returned$Days.Shelter)



#work with ppp
library(spatstat)
library(tidyverse)
#load data
setwd("/home/lionel/Desktop/Blog/Stage_2017_terec/")
dat <- read.table("data/Dataset mieren af.csv",sep=";",head=TRUE)
#coords in meters
dat$X <- dat$X/100
dat$Y <- dat$Y/100

dat <- dat[!duplicated(dat),]

#the ants ppp
ant_pp <- ppp(subset(dat,Soort=="Tetramorium_caespitum")[,"X"],subset(dat,Soort=="Tetramorium_caespitum")[,"Y"],owin(c(0,15),c(0,15)))

#put all plant into a list of ppp
lpl_pp <- lapply(levels(dat$Soort),function(sp){
  tmp <- subset(dat,Soort==sp)
  return(ppp(tmp$X,tmp$Y,owin(c(0,15),c(0,15))))
})

#turn this into a list of dens
ldens <- lapply(lpl_pp,function(x) density(x))
names(ldens) <- levels(dat$Soort)
#remove the ants
ldens <- ldens[-15]

#a density plot for all plant taxa
par(mfrow=c(4,4),mar=c(0,0,3,0))
plot(density(ant_pp),main="Ants")
for(i in 1:15){
  plot(ldens[[i]],main=names(ldens)[i])
}

#a first model with a few plant species with similar density
mm <- kppm(ant_pp,trend=~Crepis_sp.+Hypochaeris_radicata+Oenothera_sp.+Sedum_acre,data=ldens)

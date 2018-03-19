#load the functions
setwd("~/Desktop/Blog/GitFolder/AddingSE/")
lisf <- list.files(pattern="*.r",path = "functions/")
sapply(lisf,function(x) source(paste0("functions/",x)))

#code
#some simulated data
dat <- data.frame(income_std=runif(100,-2,2),area=gl(n = 4,k=25,labels = c("North","East","South","West")))
modmat <- model.matrix(~income_std*area,dat)
coeff <- c(1,2,-0.5,1.2,0.6,-1.3,-0.2,2)
dat$tree_nb <- rnorm(100,mean = modmat%*%coeff,sd=1)

#fit the model 
(m <- lm(tree_nb~income_std*area,dat))

#compute the summed SE
sum_se <- rbind(add_se(m,name_f = "area"),add_se(m,name_f = "area",name_x = "income_std"))

#a plot
library(viridis)
col_vec <- viridis(4)
plot(tree_nb~income_std,dat,pch=16,col=col_vec[dat$area])
abline(a = coef(m)[1],b=coef(m)[2],col=col_vec[1],lwd=2)
abline(a = sum_se[1,1],b=sum_se[4,1],col=col_vec[2],lwd=2)
abline(a = sum_se[2,1],b=sum_se[5,1],col=col_vec[3],lwd=2)
abline(a = sum_se[3,1],b=sum_se[6,1],col=col_vec[4],lwd=2)

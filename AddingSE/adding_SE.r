#adding standard error blog post

#function parameter
#@model: a lm or glm fitted model
#@name_f: a character string with the name of the categorical (factor) variable
#@name_x: a character string with the name fo the interacting variable, by default the intercept
add_se <- function(model,name_f,name_x="Intercept"){
  #grab the standard error of the coefficients
  se_vec <- summary(model)$coefficients[,2]
  if(name_x=="Intercept"){
    #the stabdard error of the intercept
    se_x <- se_vec[1]
    #get the level-specific standard errors
    se_f <- se_vec[grep(name_f,names(se_vec))]
    se_f <- se_f[-grep(":",names(se_f))]
    #get the covariance between the intercept and the level-specific parameters
    vcov_f <- vcov(model)[grep(name_f,rownames(vcov(model))),grep(name_x,colnames(vcov(model)))]
    vcov_f <- vcov_f[-grep(":",names(vcov_f))]
    #the estimated average value at each level
    coef_f <- coef(model)[1]+coef(model)[names(vcov_f)]
  }
  else{
    #similar code for the case of another variable than the intercept
    se_x <- se_vec[name_x]
    se_f <- se_vec[grep(name_f,names(se_vec))]
    se_f <- se_f[grep(":",names(se_f))]
    vcov_f <- vcov(model)[grep(name_f,rownames(vcov(model))),grep(name_x,colnames(vcov(model)))][,1]
    vcov_f <- vcov_f[grep(":",names(vcov_f))]
    coef_f <- coef(model)[name_x]+coef(model)[names(vcov_f)]
  }
  #compute the summed SE
  se_f <- sqrt(se_x**2+se_f**2+2*vcov_f)
  #create the output dataframe
  out <- data.frame(Coef=coef_f,SE=se_f)
  return(out)
}


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

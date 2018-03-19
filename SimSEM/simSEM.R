library(piecewiseSEM)
library(lavaan)
library(plyr)
library(dplyr)
library(ggplot2)
library(igraph)
library(dagitty)
library(AICcmodavg)


path_to_function <- "~/Desktop/Blog/GitFolder/SimSEM/"
setwd(path_to_function)

source("simSEM_function.R")
#test over different data generation type (Type), sample size (N), number of variables (X),
#varying signal strength (sd_eff) and noise importance (sd_red).
#This will return the full graphs shown in the Appendix of the manuscript.
#main figures use sd_eff fixed at 2.5 and sd_res fixed at 1. 

#the different parameter set
sims <- expand.grid(Type=c("random","exact","shuffled","complex","simple"),N=seq(20,100,20),X=c(5,7,10),C=0.3,lv=c(FALSE,TRUE),sd_eff = c(1, 2.5, 5),sd_res = c(0.5, 1, 2))

#for testing purposes
#sims <- sims[sample(1:1170,10,replace = FALSE),]
#for the random data type, no need to vary sd_eff idependenlty
sims <- sims[!(sims$Type == "random" & sims$sd_eff %in% c(2.5,5)),]
#run the code, see the code in the simSEM_function 
res_pre <- sim_for(sims)
#aggregate over the replications per parameter set
res_pre %>%
  group_by(pkg,type,Exp,N,X,Sd_eff,Sd_res)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=median(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) -> res_all

#save the file
write.table(res_all,"simSem_res.csv",row.names=FALSE)


#res_all <- read.csv("~/Desktop/Blog/GitFolder/SimSEM/simSem_res.csv",sep=" ")
res_all$Signal <- with(res_all,paste("sd_eff:",Sd_eff,"\nsd_res:",Sd_res,sep=" "))


#make the figures
gg1 <- ggplot(res_all,aes(x=N,y=PropSigM,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  labs(x="Sample size",y="Proportion of accepted models (p > 0.05)",
       title="Proportion of accepted model for different data generation")

 ggsave("prop_acc_A.eps",gg1,width=30,height=15,unit="cm",dpi=100)
 
 gg2 <- ggplot(res_all,aes(x=N,y=PropSigPaths,color=type,linetype=pkg))+geom_path()+
   facet_grid(X~Signal) +
   labs(x="Sample size",y="Proportion of significant paths (p < 0.05)",
        title="Proportion of significant paths for different data generation")
 
 ggsave("prop_path_A.eps",gg2,width=30,height=15,unit="cm",dpi=100)
# 
 gg3 <- ggplot(subset(res_all,type!="random"),aes(x=N,y=AvgRSq,color=type,linetype=pkg))+geom_path()+
   facet_grid(X~Signal) +
   labs(x="Sample size",y="Average R-square of the constituing models",
        title="Average R-square for different data generation")
# 
 ggsave("avg_rsq_A.eps",gg3,width=30,height=15,unit="cm",dpi=100)
# 
 gg4 <- ggplot(res_all,aes(x=N,y=PropCond,color=type,linetype=pkg))+geom_path()+
   facet_grid(X~Signal) +
   labs(x="Sample size",y="Average proportion of breached conditional independence",
        title="Average failed conditional independence for different data generation")
 
 ggsave("prop_cond_A.eps",gg4,width=30,height=15,unit="cm",dpi=100)

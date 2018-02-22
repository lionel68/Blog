library(piecewiseSEM)
library(lavaan)
library(plyr)
library(dplyr)
library(ggplot2)
library(igraph)
library(dagitty)


path_to_function <- "~/Desktop/Blog/GitFolder/SimSEM/"
setwd(path_to_function)

source("simSEM_function.R")
#test over different sample size (N) and number of variables (X) with random data generation
#ie there is no simulated signal in the data
sims <- expand.grid(N=seq(20,100,10),X=5:10,C=0.3,lv=c(FALSE,TRUE))

res_random <- sim_for(sims,"random")

res_random%>%
  group_by(pkg,Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=mean(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) ->res_random_dd
#write.table(res_dd,"res_dd.csv",row.names=FALSE)

#now simulate the exact distribution for the fitted model
res_exact <- sim_for(sims,"rexact")

res_exact%>%
  group_by(pkg,Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=mean(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) ->res_exact_dd
#write.table(res_dd,"res_dd.csv",row.names=FALSE)

#now for the shuffled scenario
res_shuffled <- sim_for(sims,"shuffled")

res_shuffled%>%
  group_by(pkg,Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=mean(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) ->res_shuffled_dd
#write.table(res_dd,"res_dd.csv",row.names=FALSE)

#now for the simple scenario
res_simple <- sim_for(sims,"simple")

res_simple%>%
  group_by(pkg,Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=mean(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) ->res_simple_dd
#write.table(res_dd,"res_dd.csv",row.names=FALSE)

#now for the complex scenario
res_complex <- sim_for(sims,"complex")

res_complex%>%
  group_by(pkg,Exp,N,X,type)%>%
  summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=mean(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  mutate(PropSigPaths=NbSigPaths/NbPaths) ->res_complex_dd
#write.table(res_dd,"res_dd.csv",row.names=FALSE)

#make one graph for all
res_all <- rbind(res_random_dd,res_exact_dd,res_shuffled_dd,res_complex_dd,res_simple_dd)
res_all$Data <- rep(c("Random","Exact","Shuffled","Complex","Simple"),each=108)

#res_all$type[res_all$type=="lv"] <- "lavaan"

write.table(res_all,"simSem_res.csv",row.names=FALSE)




gg1 <- ggplot(res_all,aes(x=N,y=PropSigM,color=Data,linetype=type))+geom_path()+
  facet_wrap(~X,labeller = label_both) +
  labs(x="Sample size",y="Proportion of accepted models (p > 0.05)",
       title="Proportion of accepted model for different data generation")

ggsave("prop_acc.eps",gg1,width=20,height=15,unit="cm",dpi=100)

gg2 <- ggplot(res_all,aes(x=N,y=PropSigPaths,color=Data,linetype=type))+geom_path()+
  facet_wrap(~X,labeller = label_both) +
  labs(x="Sample size",y="Proportion of significant paths (p < 0.05)",
       title="Proportion of significant paths for different data generation")

ggsave("prop_path.eps",gg2,width=20,height=15,unit="cm",dpi=100)

gg3 <- ggplot(subset(res_all,Data!="Random"),aes(x=N,y=AvgRSq,color=Data,linetype=type))+geom_path()+
  facet_wrap(~X,labeller = label_both) +
  labs(x="Sample size",y="Average R-square of the constituing models",
       title="Average R-square for different data generation")

ggsave("avg_rsq.eps",gg3,width=20,height=15,unit="cm",dpi=100)

gg4 <- ggplot(res_all,aes(x=N,y=PropCond,color=Data,linetype=type))+geom_path()+
  facet_wrap(~X,labeller = label_both) +
  labs(x="Sample size",y="Average proportion of breached conditional independence",
       title="Average failed conditional independence for different data generation")

ggsave("prop_cond.eps",gg4,width=20,height=15,unit="cm",dpi=100)
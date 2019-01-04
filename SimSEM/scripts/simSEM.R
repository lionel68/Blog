##################

# R script to reproduce the results in the submitted manuscript:

# How robust are Structural Equation Models to model mis-specification? A simulation study

# author: Lionel Hertzog

# date: 17.10.2018

########

# load the libraries
library(piecewiseSEM) # should be a version below 2.0 !!!!!!
library(lavaan)
library(igraph)
library(dagitty)
library(AICcmodavg)
library(plyr)
library(tidyverse)



# set the wd and load the helper functions
path_to_function <- "~/Desktop/Blog/GitFolder/SimSEM/"
setwd(path_to_function)

source("scripts/simSEM_function.R")

# create the parameter set to run the simulations on
# test over different data generation type (Type), sample size (N), number of variables (X),
# varying signal strength (sd_eff) and noise importance (sd_red).

# will also need to re-run this with new generation of coefficients

sims <- expand.grid(type=c("random","exact","shuffled","overspecified","underspecified"),
                    N=c(seq(20,100,20),200,500,1000,5000,10000),
                    X=c(5,7,10),
                    C=0.3,
                    p_var = 0.25,
                    pkg=c("lavaan","pcsem"),
                    sd_eff = c(1, 2.5, 5),
                    sd_res = c(0.5, 1, 2),stringsAsFactors = FALSE)


# for the random data type, no need to vary sd_eff idependenlty
sims <- sims[!(sims$type == "random" & sims$sd_eff %in% c(2.5,5)),]

# run the simulations, see the code in the simSEM_function
res_pre <- sim_for(sims,n_rep=100)

# write the output
write.table(res_pre,"data/simSem_res_pre_rev.csv",row.names=FALSE)

# load back the simulation data
res_pre <- read.csv("data/simSem_res_pre_rev.csv", sep = " ")

# aggregate over the replications per parameter set
res_pre %>%
  group_by(pkg,type,Exp,N,X,Sd_eff,Sd_res)%>%
  dplyr::summarise(PropSigM=sum(sigM,na.rm=TRUE)/n(),NbPaths = sum(nbPaths),NbSigPaths=sum(nbSigPaths),
            AvgRSq=median(avgRsq,na.rm=TRUE),PropCond = mean(Local,na.rm=TRUE))%>%
  dplyr::mutate(PropSigPaths=NbSigPaths/NbPaths) %>%
  dplyr::mutate(Signal = paste0("sd_eff: ",Sd_eff,"\nsd_res: ",Sd_res)) -> res_all

# the main figures:

## fig.1 variation in proportion of accepted models
gg1 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=PropSigM,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of accepted models (p > 0.05)",
       title="") +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/prop_acc_rev.png",gg1,width=30,height=15,unit="cm",dpi=100)
# discuss better detecting power for pcsem as model get more complex for low sample size

## fig.2 variation in proportion of accepted paths
gg2 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=PropSigPaths,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  geom_hline(yintercept = 0.8,linetype="dashed") +
  labs(x="Sample size",y="Proportion of significant paths (p < 0.05)",
       title="") +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/prop_path_rev.png",gg2,width=30,height=15,unit="cm",dpi=100)
# discuss the heuristic of good power to detect sig paths for 10 data points per paths

## fig.3 variation in average R-square
gg3 <- ggplot(subset(res_all,Sd_res == 1 & Sd_eff == 2.5),aes(x=N,y=AvgRSq,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Average R-square of the constituing models",
       title="") +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/avg_rsq_rev.png",gg3,width=30,height=15,unit="cm",dpi=100)
#
## fig.4 variation in proportion of failed conditional independence tests
gg4 <- ggplot(subset(res_all,(Sd_res == 1 & Sd_eff == 2.5) | (type == "random" & Sd_res == 1)),aes(x=N,y=PropCond,color=type,linetype=pkg))+geom_path()+
  facet_grid(~X) +
  scale_x_log10() +
  labs(x="Sample size",y="Average proportion of breached conditional independence",
       title="") +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/prop_cond_rev.png",gg4,width=30,height=15,unit="cm",dpi=100)

#
# # appendix figures
# ## this time looking across signal / noise ratios
#
gg1 <- ggplot(res_all,aes(x=N,y=PropSigM,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Proportion of accepted models (p > 0.05)",
       title="Proportion of accepted model for different data generation") +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/prop_acc_A_rev.png",gg1,width=30,height=15,unit="cm",dpi=100)

gg2 <- ggplot(res_all,aes(x=N,y=PropSigPaths,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  geom_hline(yintercept = 0.8,linetype="dashed") +
  labs(x="Sample size",y="Proportion of significant paths (p < 0.05)",
       title="Proportion of significant paths for different data generation") +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/prop_path_A_rev.png",gg2,width=30,height=15,unit="cm",dpi=100)
# #
gg3 <- ggplot(subset(res_all,type!="random"),aes(x=N,y=AvgRSq,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Average R-square of the constituing models",
       title="Average R-square for different data generation")  +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")
#
ggsave("figures/avg_rsq_A_rev.png",gg3,width=30,height=15,unit="cm",dpi=100)
# #
gg4 <- ggplot(res_all,aes(x=N,y=PropCond,color=type,linetype=pkg))+geom_path()+
  facet_grid(X~Signal) +
  scale_x_log10() +
  labs(x="Sample size",y="Average proportion of breached conditional independence",
       title="Average failed conditional independence for different data generation")  +
  scale_linetype_discrete(name = "package") +
  scale_color_discrete(name= "data generation\nscenario")

ggsave("figures/prop_cond_A_rev.png",gg4,width=30,height=15,unit="cm",dpi=100)
#


# second set of simulations to compute and compare IC

# new set of sims with more accurate IC comparisons
sims_c <- expand.grid(N = c(seq(20,100,20),200,500,1000,5000,10000),
                    X = c(5,7,10),
                    C = 0.3,
                    p_var = 0.25,
                    pkg = c("lavaan","pcsem"),
                    sd_eff = c(1, 2.5, 5),
                    sd_res = c(0.5, 1, 2))
# 
# 
# # run the simulations, see the code in the simSEM_function 
# n_rep <- 100
# res_pre <- sim_for_c(sims_c,n_rep=n_rep) # still some bug lurking about all endogenous variables being conditionaly dependent
# 
# write.table(res_pre,"data/res_pre_c2.csv",row.names=FALSE)
# 
# # next step would be to get how often each scenarios are identified as the best
# # so comparing for each rows the best model leading value such as: unspecified, exact, shuffled or overspecfied, then aggregate this per Exp as prop
# 
# load back saved data
res_pre <- read.csv("data/res_pre_c2.csv",sep=" ")

sims_c$Exp <- 1:nrow(sims_c)

res_pre %>%
  mutate(rep = 1:n()) %>% # add an unique ID for each replicate
  gather("Type","value",AIC_t:hbic_o) %>% # melt
  separate(Type,c("metric","scenario")) %>% # separate IC metric from the scenario
  group_by(rep,N,X,sd_eff,sd_res,C,p_var,metric) %>% # group by replicate ID
  mutate(min = min(value)) %>% # the minimum value per replicate and metric
  mutate(is_small = if_else(abs(value - min) < 2,1,0)) %>% # if the IC value is within 2 unit of the smallest one mark as 1
  summarise(best_scenario = if_else(sum(is_small) == 1,scenario[which.min(value)],"none")) %>% # if there is only one smallest IC value get the scenario
  ungroup() %>%
  mutate(Exp = rep(1:nrow(sims_c),each = 3 * n_rep)) %>% # re-add the Exp ID
  group_by(Exp, N, X, sd_eff, sd_res, C, p_var, metric) %>% # group by Exp ID (ie rows in the sims_c dataframe)
  count(best_scenario) %>% # count how often each scenario pops out across the replicate within Exp ID
  mutate(prop = n / sum(n)) %>%
  left_join(sims_c[,c("Exp","pkg")],by="Exp") -> dd # get the proportion

dd$best_scenario <- factor(dd$best_scenario,labels = c("none","overspecified","shuffled","exact","underspecified"))

dd_aic <- subset(dd,metric=="AIC" & sd_eff == 2.5 & sd_res == 1)

gg1 <- ggplot(dd_aic,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(~X) +
  scale_x_log10() +
  labs(title="AICc")


dd_bic <- subset(dd,metric=="BIC" & sd_eff == 2.5 & sd_res == 1)

gg2 <- ggplot(dd_bic,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(~X) +
  scale_x_log10() +
  labs(title="BIC")

dd_hbic <- subset(dd,metric=="hbic" & sd_eff == 2.5 & sd_res == 1)

gg3 <- ggplot(dd_hbic,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(~X) +
  scale_x_log10() +
  labs(title="hbic")

dd_all <- subset(dd, sd_eff ==2.5 & sd_res == 1 & !is.na(dd$best_scenario))

gg_all <- ggplot(dd_all,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(X~metric) +
  scale_x_log10(breaks=c(20,50,100,1000,10000)) +
  labs(title="",x="Sample size",y="Proportion of simulations") +
  scale_linetype_discrete(name="package") +
  scale_color_discrete(name="best scenario") +
  theme(panel.grid.minor = element_blank())

ggsave("figures/ic_metric_prop.png",gg_all,width=10,height=8)

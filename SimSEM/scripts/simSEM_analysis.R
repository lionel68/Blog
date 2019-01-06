##################

# R script to plot the results from the simulations from the submitted manuscript:

# How robust are Structural Equation Models to model mis-specification? A simulation study

# author: Lionel Hertzog

# date: 06.01.2019

########

# note that you can re-run the simulations on your machine using the scripts simSEM_simulations.R and simSEM_function.R, you can also download the data used in the manuscript from: URL link

## load the libraries
library(tidyverse)

## set the wd, point to where the github repo was cloned
setwd("")

## first simulation batch

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

## second simulation batch

# load back saved data
res_pre <- read.csv("data/res_pre_c2.csv",sep=" ")

# the parameter set used to generate the simulations
sims_c <- expand.grid(N = c(seq(20,100,20),200,500,1000,5000,10000),
                    X = c(5,7,10),
                    C = 0.3,
                    p_var = 0.25,
                    pkg = c("lavaan","pcsem"),
                    sd_eff = c(1, 2.5, 5),
                    sd_res = c(0.5, 1, 2))
# add an ID to each line
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
  left_join(sims_c[,c("Exp","pkg")],by="Exp") %>% # get the proportion
  mutate(Signal = paste0("sd_eff: ",Sd_eff,"\nsd_res: ",Sd_res)) -> dd 

dd$best_scenario <- factor(dd$best_scenario,labels = c("none","overspecified","shuffled","exact","underspecified"))

## figure 5. proportion of replicates where each scenario was labelled best scenario
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

## appendix figures looking across signal / noise ratio for the different IC metrics separately
dd_aic <- subset(dd,metric=="AIC")

gg_aic <- ggplot(dd_aic,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(Signal~X) +
  scale_x_log10() +
  labs(title="AICc")

ggsave("figures/aic_appendix.png",gg_aic, width = 10, height = 8)


dd_bic <- subset(dd,metric=="BIC" & sd_eff == 2.5 & sd_res == 1)

gg_bic <- ggplot(dd_bic,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(Signal~X) +
  scale_x_log10() +
  labs(title="BIC")

ggsave("figures/bic_appendix.png",gg_bic,height=8,width=10)

dd_hbic <- subset(dd,metric=="hbic" & sd_eff == 2.5 & sd_res == 1)

gg_hibc <- ggplot(dd_hbic,aes(x=N,y=prop,color=best_scenario,linetype=pkg)) +
  geom_line() +
  facet_grid(Signal~X) +
  scale_x_log10() +
  labs(title="hbic")

ggsave("figures/hbic_appendix.png",gg_hbic, width = 10, height = 8)

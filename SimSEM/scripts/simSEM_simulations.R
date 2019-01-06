##################

# R script to run the simulations from the submitted manuscript:

# How robust are Structural Equation Models to model mis-specification? A simulation study

# author: Lionel Hertzog

# date: 06.01.2019

########

# load the libraries
library(piecewiseSEM) # should be a version below 2.0 !!!!!!
library(lavaan)
library(igraph)
library(dagitty)
library(AICcmodavg)
library(plyr)



# set the wd, point to where the github repo was cloned
path_to_function <- "~/Desktop/Blog/GitFolder/SimSEM/"
setwd(path_to_function)

# load the functions
source("scripts/simSEM_function.R")

## first simulation batch

# create the parameter set to run the simulations on
# test over different data generation type (type), sample size (N), number of variables (X),
# varying signal strength (sd_eff), noise importance (sd_red) and estimation method (pkg).
# Constant parameters (needed in the functions) are: connectance of the DAG (C) and
# the proportion of paths to vary in data generation shuffled, overspecified or
# underspecified (p_var)

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

# run the simulations, with 100 replication per parameter set
# see the code in the simSEM_function for details
res_pre <- sim_for(sims,n_rep=100)

# write the output
write.table(res_pre,"data/simSem_res_pre_rev.csv",row.names=FALSE)


## second simulation batch focusing on Information Criteria

# create the parameter set
sims_c <- expand.grid(N = c(seq(20,100,20),200,500,1000,5000,10000),
                    X = c(5,7,10),
                    C = 0.3,
                    p_var = 0.25,
                    pkg = c("lavaan","pcsem"),
                    sd_eff = c(1, 2.5, 5),
                    sd_res = c(0.5, 1, 2))
 
# run the simulations, see the code in the simSEM_function 
n_rep <- 100
res_pre <- sim_for_c(sims_c,n_rep=n_rep) 
 
 write.table(res_pre,"data/res_pre_c2.csv",row.names=FALSE)



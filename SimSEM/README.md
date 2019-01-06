# Code for: How robusts are Structural Equation Models to model miss-specification

An arXiv version of the manuscript can be found [here](https://arxiv.org/abs/1803.06186). This manuscript explore how different model miss-specification affect various SEM fitness metrics.

In this repo is stored the code used to generate the figures shown in the manuscript. There are three R files to fully reproduce these results in the scripts folder:

* simSEM_simulations.R: contain the code to run the simulations

* simSEM_function.R: contain all the functions used to run the simulations

* simSEM_analysis.R: contain the code to produce the figures shown in the manuscript

Alternatively the data used to produce the figures shown in the manuscript can be downloaded from: 

The parameters that can be varied in the simulations are:

* Type: the type of data-generation, one character from the following, random, exact, overspecified, underspecified or shuffled.
* X: the number of covariates in the models
* N: the sample size
* C: the graph connectance, how many links are drawn between the variables
* pkg: whether to fit the model via lavaan (lavaan) or piecewiseSEM (pcsem)
* sd_res: the standard deviation of the residuals, ie the amount of white noise
* sd_eff: the standard deviation of the path coefficients, ie the amount of signal, larger values means that larger path coefficients are more likely to be generated.

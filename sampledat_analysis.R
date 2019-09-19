source(file="~/codes_submission/fun_general.R")

library(stats)
library(MASS)
library(mvtnorm)
library(mice)
library(glmnet)
require(e1071)
library(ncvreg)

## READ IN TRAINING/TEST DATA

data.obs<-read.csv("~/codes_submission/traindat.csv")

#####################################
#### Var Sel: BISS   ################
#####################################

alph<-0.5
reps<-800
reps2<-500
thresh<-seq(1,0.1,-0.025)
lambda.prop=0.75
m=25 # we have 25 markers
K=1 # we have only one W


set.seed(1.234)

biss_p_log<-BI_SS(data.obs, impute_freq_gauss, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop)
beta.biss.p.las<-biss_p_log$beta.las
beta.biss.p.las

biss_mice_log<-BI_SS(data.obs, impute_freq_mice, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop)
beta.biss.m.las<-biss_mice_log$beta.las
beta.biss.m.las

biss_n_log<-BI_SS(data.obs, impute_freq_nonpara, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop)
beta.biss.n.las<-biss_n_log$beta.las
beta.biss.n.las

biss_e_log<-BI_SS(data.obs, impute_freq_exp, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop)
beta.biss.e.las<-biss_e_log$beta.las
beta.biss.e.las

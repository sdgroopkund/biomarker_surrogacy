
##############################################
#### Variable Selection steps          #######
##############################################

## NOTE: Need to simulate data before running this script (see simu codes).

## NOTE that the variable weights, if present in a dataset, should be the sampling weights. 
## The methods below implement the inverse of these sampling weights for the weighted analyses.
## Current functions can only treat the case when sampling weights vary only with Y and Z


if (!is.null(data.obs$weights)) data.obs$weights<-rep(1,nrow(data.obs)) # sampling weights

#######################################
#### Var Sel 1: no selection  #########
#######################################

n.multimp<-10

snam <- paste(paste0("I(W*S.", 1:m,")",collapse="+"),paste0("I(Z*W*S.",1:m,")",collapse="+"),sep='+')
fmla <- as.formula(paste("Y~Z+",snam,"+W+I(W*Z)"))
beta.logit_p<-matrix(0,n.multimp,2*(m+K+1))
beta.logit_m<-matrix(0,n.multimp,2*(m+K+1))
beta.logit_n<-matrix(0,n.multimp,2*(m+K+1))
beta.logit_e<-matrix(0,n.multimp,2*(m+K+1))
for (k in 1:n.multimp){
  data.imp_p<-impute_freq_gauss(data.obs,m)$data.imp
  fit.logit<-glm(fmla,data=data.imp_p,weights = 1/data.imp_p$weights,family=binomial(link=logit)) ## RISK MODEL LOGIT
  beta.logit_p[k,]<-fit.logit$coef
  
  data.imp_mice<-impute_freq_mice(data.obs,m)$data.imp
  fit.logit<-glm(fmla,data=data.imp_mice,weights = 1/data.imp_mice$weights,family=binomial(link=logit)) 
  beta.logit_m[k,]<-fit.logit$coef
  
  data.imp_n<-impute_freq_nonpara(data.obs,m)$data.imp
  fit.logit<-glm(fmla,data=data.imp_n,weights = 1/data.imp_n$weights,family=binomial(link=logit)) 
  beta.logit_n[k,]<-fit.logit$coef
  
  data.imp_e<-impute_freq_exp(data.obs,m)$data.imp
  fit.logit<-glm(fmla,data=data.imp_e,weights = 1/data.imp_e$weights,family=binomial(link=logit)) 
  beta.logit_e[k,]<-fit.logit$coef
  
}
beta.nosel.p<-colMeans(beta.logit_p)
beta.nosel.m<-colMeans(beta.logit_m)
beta.nosel.n<-colMeans(beta.logit_n)
beta.nosel.e<-colMeans(beta.logit_e)

#######################################
#### Var Sel 2: DIRECT ################
#######################################



beta.logit_p_las<-matrix(0,n.multimp,2*(m+K+1))
beta.logit_m_las<-matrix(0,n.multimp,2*(m+K+1))
beta.logit_n_las<-matrix(0,n.multimp,2*(m+K+1))
beta.logit_e_las<-matrix(0,n.multimp,2*(m+K+1))
for (k in 1:n.multimp){
  # PARA GAUSS
  data.imp_p<-impute_freq_gauss(data.obs,m)$data.imp
  indices<-which(apply(data.imp_p,1,function(x) sum(is.na(x)))>0)
  data.imp_p<-data.imp_p[setdiff(1:nrow(data.imp_p),indices),]
  datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                          paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp_p)
  datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
  cv.las<-cv.glmnet(as.matrix(datamat[,-1]),data.imp_p$Y,weights = 1/data.imp_p$weights,family='binomial',standardize=F)
  fit.las<-glmnet(as.matrix(datamat[,-1]),data.imp_p$Y,weights = 1/data.imp_p$weights,lambda=cv.las$lambda.min,family='binomial',standardize=F)
  beta.las_p<-as.numeric(fit.las$beta)
  datamat.las<-as.data.frame(datamat[,which(beta.las_p!=0),drop=F])
  if(1%in%which(beta.las_p!=0)) datamat.las.g=datamat.las[,-1,drop=F] else datamat.las.g=datamat.las
  datamat.las.g$Y<-data.imp_p$Y
  beta.logit_p_las[k,union(c(1),which(beta.las_p!=0))]<-glm(Y~.,data=datamat.las.g,weights = 1/data.imp_p$weights,family=binomial(link='logit'))$coef

  #MICE
  data.imp_m<-impute_freq_mice(data.obs,m)$data.imp
  indices<-which(apply(data.imp_m,1,function(x) sum(is.na(x)))>0)
  data.imp_m<-data.imp_m[setdiff(1:nrow(data.imp_m),indices),]
  datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                          paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp_m)
  datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
  cv.las<-cv.glmnet(as.matrix(datamat[,-1]),data.imp_m$Y,weights = 1/data.imp_m$weights,family='binomial',standardize=F)
  fit.las<-glmnet(as.matrix(datamat[,-1]),data.imp_m$Y,weights = 1/data.imp_m$weights,lambda=cv.las$lambda.min,family='binomial',standardize=F)
  beta.las_m<-as.numeric(fit.las$beta)
  datamat.las<-as.data.frame(datamat[,which(beta.las_m!=0),drop=F])
  if(1%in%which(beta.las_m!=0)) datamat.las.g=datamat.las[,-1,drop=F] else datamat.las.g=datamat.las
  datamat.las.g$Y<-data.imp_m$Y
  beta.logit_m_las[k,union(c(1),which(beta.las_m!=0))]<-glm(Y~.,data=datamat.las.g,weights = 1/data.imp_m$weights,family=binomial(link='logit'))$coef
  
  #NONPARA
  data.imp_n<-impute_freq_nonpara(data.obs,m)$data.imp
  indices<-which(apply(data.imp_n,1,function(x) sum(is.na(x)))>0)
  data.imp_n<-data.imp_n[setdiff(1:nrow(data.imp_n),indices),]
  datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                          paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp_n)
  datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
  cv.las<-cv.glmnet(as.matrix(datamat[,-1]),data.imp_n$Y,weights = 1/data.imp_n$weights,family='binomial',standardize=F)
  fit.las<-glmnet(as.matrix(datamat[,-1]),data.imp_n$Y,weights = 1/data.imp_n$weights,lambda=cv.las$lambda.min,family='binomial',standardize=F)
  beta.las_n<-as.numeric(fit.las$beta)
  datamat.las<-as.data.frame(datamat[,which(beta.las_n!=0),drop=F])
  if(1%in%which(beta.las_n!=0)) datamat.las.g=datamat.las[,-1,drop=F] else datamat.las.g=datamat.las
  datamat.las.g$Y<-data.imp_n$Y
  beta.logit_n_las[k,union(c(1),which(beta.las_n!=0))]<-glm(Y~., data=datamat.las.g, weights = 1/data.imp_n$weights,family=binomial(link='logit'))$coef
  
  #PARA EXP
  data.imp_e<-impute_freq_exp(data.obs,m)$data.imp
  indices<-which(apply(data.imp_e,1,function(x) sum(is.na(x)))>0)
  data.imp_e<-data.imp_e[setdiff(1:nrow(data.imp_e),indices),]
  datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                          paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp_e)
  datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
  cv.las<-cv.glmnet(as.matrix(datamat[,-1]),data.imp_e$Y,weights = 1/data.imp_e$weights,family='binomial',standardize=F)
  fit.las<-glmnet(as.matrix(datamat[,-1]),data.imp_e$Y,weights = 1/data.imp_e$weights,lambda=cv.las$lambda.min,family='binomial',standardize=F)
  beta.las_e<-as.numeric(fit.las$beta)
  datamat.las<-as.data.frame(datamat[,which(beta.las_e!=0),drop=F])
  if(1%in%which(beta.las_e!=0)) datamat.las.g=datamat.las[,-1,drop=F] else datamat.las.g=datamat.las
  datamat.las.g$Y<-data.imp_e$Y
  beta.logit_e_las[k,union(c(1),which(beta.las_e!=0))]<-glm(Y~.,data=datamat.las.g,weights = 1/data.imp_e$weights,family=binomial(link='logit'))$coef
}


beta.direct.p.las<-colMeans(beta.logit_p_las)
beta.direct.m.las<-colMeans(beta.logit_m_las)
beta.direct.n.las<-colMeans(beta.logit_n_las)
beta.direct.e.las<-colMeans(beta.logit_e_las)

#######################################
#### Var Sel 3: BISS   ################
#######################################

alph<-0.5
reps<-800
reps2<-500
thresh<-seq(1,0.1,-0.025)
lambda.prop=0.75

biss_p_log<-BI_SS(data.obs, impute_freq_gauss, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop = lambda.prop)
beta.biss.p.las<-biss_p_log$beta.las

biss_mice_log<-BI_SS(data.obs, impute_freq_mice, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop = lambda.prop)
beta.biss.m.las<-biss_mice_log$beta.las

biss_n_log<-BI_SS(data.obs, impute_freq_nonpara, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop = lambda.prop)
beta.biss.n.las<-biss_n_log$beta.las

biss_e_log<-BI_SS(data.obs, impute_freq_exp, alph, link='logit', reps, reps2, m, K, thresh, lambda.prop = lambda.prop)
beta.biss.e.las<-biss_e_log$beta.las


##############################################
####    Predicting with each beta       ######
##############################################

pred.nosel.p<-predict.mod(datamat.t,datain.t$Y,beta.nosel.p,'logit','logit',sig,beta)$scores
pred.nosel.n<-predict.mod(datamat.t,datain.t$Y,beta.nosel.n,'logit','logit',sig,beta)$scores
pred.nosel.e<-predict.mod(datamat.t,datain.t$Y,beta.nosel.e,'logit','logit',sig,beta)$scores
pred.nosel.m<-predict.mod(datamat.t,datain.t$Y,beta.nosel.m,'logit','logit',sig,beta)$scores
pred.direct.p.las<-predict.mod(datamat.t,datain.t$Y,beta.direct.p.las,'logit','logit',sig,beta)$scores
pred.direct.n.las<-predict.mod(datamat.t,datain.t$Y,beta.direct.n.las,'logit','logit',sig,beta)$scores
pred.direct.e.las<-predict.mod(datamat.t,datain.t$Y,beta.direct.e.las,'logit','logit',sig,beta)$scores
pred.direct.m.las<-predict.mod(datamat.t,datain.t$Y,beta.direct.m.las,'logit','logit',sig,beta)$scores

pred.biss.p.las<-matrix(0,8,length(thresh))
pred.biss.n.las<-matrix(0,8,length(thresh))
pred.biss.m.las<-matrix(0,8,length(thresh))
pred.biss.e.las<-matrix(0,8,length(thresh))
for (t in 1:length(thresh)){
  pred.biss.p.las[,t]<-predict.mod(datamat.t,datain.t$Y,beta.biss.p.las[,t],'logit','logit',sig,beta)$scores
  pred.biss.m.las[,t]<-predict.mod(datamat.t,datain.t$Y,beta.biss.m.las[,t],'logit','logit',sig,beta)$scores
  pred.biss.n.las[,t]<-predict.mod(datamat.t,datain.t$Y,beta.biss.n.las[,t],'logit','logit',sig,beta)$scores
  pred.biss.e.las[,t]<-predict.mod(datamat.t,datain.t$Y,beta.biss.e.las[,t],'logit','logit',sig,beta)$scores
}



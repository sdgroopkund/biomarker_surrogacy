## NOTE that the variable weights, if present in a dataset, should be the sampling weights. 
## The methods below implement the inverse of these sampling weights for the weighted analyses.
## Current functions can only treat the case when sampling weights vary only with Y and Z



# 1. Expit transformation
expit<-function(x) 1/(1+exp(-x))


#########################################
###### IMPUTATION FUNCTIONS #############
#########################################

# 2. Impute with MICE
impute_freq_mice<-function(data.obs,n.mark,m=1){
  ini <- mice(data.obs, maxit=0)
  pred<-ini$pred
  mice.obj<-mice(data.obs,m,print=F,pred=pred)
  data.imp_mice<-complete(mice.obj)
  names(data.imp_mice)[1:n.mark]<-paste0('S.',1:n.mark)
  return(list(data.imp=data.imp_mice,pred=pred))
}

# 3. Impute with Parametric (GAUSS)
impute_freq_gauss<-function(dataout,m){
  
  S<-as.matrix(subset(dataout,select=1:m))   ### S(1) for marker 1,..,m
  Z<-dataout$Z
  Y<-dataout$Y  
  W<-dataout$W
  wt<-rep(1,nrow(dataout))
  if (!is.null(dataout$weights)) wt<-1/dataout$weights #inverse of sampling weights
  
  S.z0.y0<-S[Z==0 & Y==0,,drop=F]
  wt.z0.y0<-wt[Z==0 & Y==0]
  W.z0.y0<-W[Z==0 & Y==0]
  S.z0.y1<-S[Z==0 & Y==1,,drop=F]
  wt.z0.y1<-wt[Z==0 & Y==1]
  W.z0.y1<-W[Z==0 & Y==1]
  S.z1<-S[Z==1,,drop=F]
  wt.z1<-wt[Z==1]
  W.z1<-W[Z==1]
  S.z0<-S[Z==0,,drop=F]
  wt.z0<-wt[Z==0]
  W.z0<-W[Z==0]
  
  fit.z0<-glm(Y~1+W,data=dataout[dataout$Z==0,],weights=wt.z0,family=binomial(link=logit))
  a5_hat<-fit.z0$coef[1]
  a6_hat<-fit.z0$coef[2:length(fit.z0$coef)]
  
  if (is.factor(dataout$W)){
    Wun<-levels(W)
    pi_hat<-function(w) {
      pi<-rep(0,length(w))
      for(i in 1:length(w)){
        indw<-which(Wun==w[i])
        vec<-c(1,rep(0,length(Wun)-1))
        vec[indw]<-1
        pi[i]<-expit(sum(c(a5_hat,a6_hat)*vec))
      }
      return(pi)
    }
    Zmat.z0.y0<-model.matrix(~W.z0.y0)
    beta.hat.z0.y0<-solve(t(Zmat.z0.y0)%*%Zmat.z0.y0)%*%t(Zmat.z0.y0)%*%S.z0.y0
    cov.hat.z0.y0<-t(S.z0.y0-Zmat.z0.y0%*%beta.hat.z0.y0)%*%(S.z0.y0-Zmat.z0.y0%*%beta.hat.z0.y0)/(nrow(S.z0.y0)-ncol(Zmat.z0.y0)-1)
    a3_hat<-beta.hat.z0.y0[1,]
    a2_hat<-beta.hat.z0.y0[2:nrow(beta.hat.z0.y0),]
    Sigma_hat<-cov.hat.z0.y0
    Wmat.z0<-model.matrix(~W.z0)  
    
    E.S.z1.mat<-sapply(1:length(Wun),function(i) apply(S.z1[W.z1==Wun[i],],2,
                                                       function(x) wtd.mean(x,wt.z1[W.z1==Wun[i]])))
    a1_ind<- function(w){
      indw=which(Wun==w)
      mean1=E.S.z1.mat[,indw]
      pi=pi_hat(w)
      if (indw==1) mean2=0 else mean2=a2_hat[indw-1,]
      a1=(mean1-(1-pi)*a3_hat-mean2)/pi
      return(a1)
    }
    a1_mat<-matrix(0,length(W.z0),ncol(S))
    for (i in 1:length(W.z0)) a1_mat[i,]<-a1_ind(W.z0[i])
    a1_hat<-apply(a1_mat,2,function(x) wtd.mean(x,wt.z0))
    
    mu_S_z0_y1<-sapply(1:length(W[Z==0 & Y==1]), function(i) {
      indi=which(Wun==W[Z==0 & Y==1][i])
      if (indi ==1) mean=a1_hat else mean=a1_hat+a2_hat[indi-1,]
      return(mean)
    })
    
  } else{
    pi_hat<-function(w) expit(a5_hat+a6_hat%*%w)
    colm<-apply(S.z1,2,function(x) wtd.mean(x,wt.z1))
    Zmat.z0.y0<-cbind(rep(1,nrow(S.z0.y0)),W.z0.y0)
    beta.hat.z0.y0<-solve(t(Zmat.z0.y0)%*%Zmat.z0.y0)%*%t(Zmat.z0.y0)%*%S.z0.y0
    cov.hat.z0.y0<-t(S.z0.y0-Zmat.z0.y0%*%beta.hat.z0.y0)%*%(S.z0.y0-Zmat.z0.y0%*%beta.hat.z0.y0)/(nrow(S.z0.y0)-ncol(Zmat.z0.y0)-1)
    a3_hat<-beta.hat.z0.y0[1,]
    a2_hat<-beta.hat.z0.y0[2,]
    Sigma_hat<-cov.hat.z0.y0
    a1_hat<-(colm-wtd.mean(1-pi_hat(W.z0),wt.z0)*a3_hat-a2_hat*wtd.mean(W.z0,wt.z0))/wtd.mean(pi_hat(W.z0),wt.z0)
    mu_S_z0_y1<-sapply(1:length(W[Z==0 & Y==1]), function(i) a1_hat+a2_hat*W[Z==0 & Y==1][i])
    
  }
  
  S_imp<-sapply(1:length(W.z0.y1), function(i) rmvnorm(1,mu_S_z0_y1[,i],Sigma_hat))
  S[Z==0 & Y==1,]<-t(S_imp)
  data_imp<-dataout
  data_imp[,1:m]<-S
  return(list(data.imp=data_imp,a1_hat=a1_hat,a2_hat=a2_hat,a3_hat=a3_hat,a5_hat=a5_hat,Sigma_hat=Sigma_hat,pi_hat=pi_hat))
}



# 4. Impute with Parametric (Exponential) - currently runs for only continuous W
impute_freq_exp<-function(dataout,m){
  
  S<-as.matrix(subset(dataout,select=1:m))   ### S(1) for marker 1,..,m
  Z<-dataout$Z
  Y<-dataout$Y   
  W<-dataout$W    
  
  wt<-rep(1,nrow(dataout))
  if (!is.null(dataout$weights)) wt<-1/dataout$weights #inverse of sampling weights
  
  S.z0.y0<-S[Z==0 & Y==0,,drop=F]
  wt.z0.y0<-wt[Z==0 & Y==0]
  W.z0.y0<-W[Z==0 & Y==0]
  S.z0.y1<-S[Z==0 & Y==1,,drop=F]
  wt.z0.y1<-wt[Z==0 & Y==1]
  W.z0.y1<-W[Z==0 & Y==1]
  S.z1<-S[Z==1,,drop=F]
  wt.z1<-wt[Z==1]
  W.z1<-W[Z==1]
  S.z0<-S[Z==0,,drop=F]
  wt.z0<-wt[Z==0]
  W.z0<-W[Z==0]
  
  fit.z0<-glm(Y~1+W,data=dataout[dataout$Z==0,],weights=wt.z0,family=binomial(link=logit))
  a5_hat<-fit.z0$coef[1]
  a6_hat<-fit.z0$coef[2:length(fit.z0$coef)]
  pi_hat<-function(W) expit(a5_hat+a6_hat%*%W)
  
  # inverse of expectation (not weighted as sampling weights are assumed to be constant for a given level of Y and Z)
  lambda_0_inv_hat<-colMeans(S.z0.y0*W.z0.y0)
  
  library(e1071)
  library(glmnet)
  
  colm<-apply(S.z1*W.z1,2,function(x) wtd.mean(x,wt.z1))
  lambda_1_inv_hat<-pmax(colm/wtd.mean(pi_hat(W),wt)- wtd.mean(1-pi_hat(W),wt)/wtd.mean(pi_hat(W),wt)*lambda_0_inv_hat,0)
  lambda_1_hat<-1/lambda_1_inv_hat
  W.z0.y1<-W[Z==0 & Y==1]
  S_imp<-sapply(1:length(W.z0.y1), function(i){
    rexp(m, lambda_1_hat*W.z0.y1[i])
  })
  S[Z==0 & Y==1,]<-t(S_imp)
  data_imp<-dataout
  data_imp[,1:m]<-S
  return(list(data.imp=data_imp,lambda_1_hat=lambda_1_hat,lambda_0_hat=1/lambda_0_inv_hat,a5_hat=a5_hat,a6_hat=a6_hat,pi_hat=pi_hat))
}


# 5. Impute with Nonparametric 
impute_freq_nonpara<-function(dataout,m,maxit=1000){
  
  S<-as.matrix(subset(dataout,select=1:m))   ### S(1) for marker 1,..,m
  Z<-dataout$Z
  Y<-dataout$Y   
  W<-dataout$W
  
  wt<-rep(1,nrow(dataout))
  if (!is.null(dataout$weights)) wt<-1/dataout$weights #inverse of sampling weights

  if (is.factor(dataout$W)){
    Wun<-levels(W)
    pi_hat<-function(w) {
      pi<-rep(0,length(w))
      for(i in 1:length(w)){
        indw<-which(Wun==w[i])
        vec<-c(1,rep(0,length(Wun)-1))
        vec[indw]<-1
        pi[i]<-expit(sum(c(a5_hat,a6_hat)*vec))
      }
      return(pi)
    }
  } else {
    Wun<-unique(W)
    pi_hat<-function(W) {
      W<-as.matrix(W)
      return(expit(a5_hat+a6_hat%*%W))
    }
  }
  
  S.z0.y0<-S[Z==0 & Y==0,,drop=F]
  wt.z0.y0<-wt[Z==0 & Y==0]
  W.z0.y0<-W[Z==0 & Y==0]
  S.z0.y1<-S[Z==0 & Y==1,,drop=F]
  wt.z0.y1<-wt[Z==0 & Y==1]
  W.z0.y1<-W[Z==0 & Y==1]
  S.z1<-S[Z==1,,drop=F]
  wt.z1<-wt[Z==1]
  W.z1<-W[Z==1]
  S.z0<-S[Z==0,,drop=F]
  wt.z0<-wt[Z==0]
  W.z0<-W[Z==0]
  n00<-nrow(S.z0.y0)
  n1<-nrow(S.z1)
  
  fit.z0<-glm(Y~1+W,data=dataout[dataout$Z==0,],weights=wt.z0,family=binomial(link=logit))
  a5_hat<-fit.z0$coef[1]
  a6_hat<-fit.z0$coef[2:length(fit.z0$coef)]
  
  library(rootSolve)
  library(spatstat)
  
  #### Produce an estimate of the distribution functions F_1, F_00 and F_01 for a single W vector  
  
  empfun1<-vector('list',length(Wun))
  empfun00<-vector('list',length(Wun))
  for (i in 1:length(Wun)){
    S.z1.w<-S.z1[W.z1==Wun[i],]
    S.z0.y0.w<-S.z0.y0[W.z0.y0==Wun[i],]
    wt.z1.w<-wt.z1[W.z1==Wun[i]]
    wt.z0.y0.w<-wt.z0.y0[W.z0.y0==Wun[i]]
    empfun1[[i]]<-apply(S.z1.w,2,function(x) ewcdf(x,weights=wt.z1.w))
    empfun00[[i]]<-apply(S.z0.y0.w,2,function(x) ewcdf(x,weights=wt.z0.y0.w))
  }
  
  F_1_hat.w<-function(x,w){
    indw<-which(Wun==w)
    vec<-unlist(lapply(empfun1[[indw]],function(z) z(x)))
    return(vec)
  }
  
  F_00_hat.w<-function(x,w){
    indw<-which(Wun==w)
    vec<-unlist(lapply(empfun00[[indw]],function(z) z(x)))
    return(vec)
  }
  
  F_01_hat.w<-function(x,w){
    f_hat<-(F_1_hat.w(x,w)-(1-c(pi_hat(w)))*F_00_hat.w(x,w))/c(pi_hat(w))
    return(pmin(pmax(f_hat,0),1))
  }
  
  #### create a mesh of probabilities for F_01 and 
  #### Simulate missing values of S|Z=0,Y=1
  
  tstart=Sys.time()
  F_01_mesh<-vector('list',0)
  list.W<-vector('list',0)
  n.list<-0
  xmesh<-seq(range(S,na.rm=T)[1],range(S,na.rm=T)[2],length.out=10000)
  S_imp<-matrix(0,length(W.z0.y1),m)
  for (i in 1:length(W.z0.y1)){
    cont.w<-lapply(list.W,function(x) all(x==W.z0.y1[i]))
    if (sum(unlist(cont.w))==0){
      n.list<-n.list+1
      list.W[[n.list]]<-W.z0.y1[i]
      F_01_raw<-t(sapply(1:length(xmesh),function(j) F_01_hat.w(xmesh[j],W.z0.y1[i])))
      F_01_mesh[[n.list]]<-round(apply(F_01_raw,2, function(x) isoreg(xmesh,x)$yf),10)
      indi<-n.list
    } else {
      indi<-which(cont.w==T)
    }
    urand<-runif(m,0,1)
    S_imp[i,]<-sapply(1:m, function(k) {
      t<-findInterval(urand[k],sort(F_01_mesh[[indi]][,k]),all.inside=T) # this sort is not harmful because of isoreg
      xmesh[t]
    } )
   }
  tend=Sys.time()-tstart
  
  S[Z==0 & Y==1,]<-S_imp
  data_imp<-dataout
  data_imp[,1:m]<-S
  return(list(data.imp=data_imp))
}



#######################################################
###### BOOTSTRAP IMPUTATION AND STABLE SELECTION ######
#######################################################


# 6. The BISS function
BI_SS<-function(data.obs,impute.fnc,alph,link,reps,reps2,m,K,thresh,lambda.prop=NULL){
  
  # m : number of biomarkers
  # K : size of W
  # impute.fnc : the imputation function
  # alph : weight in randomized LASSO step
  # link : GLM link (probit or logit currently)
  # reps : number of bootstrap replicates in feature selection step
  # reps2 : number of bootstrap replicates in parameter estimation step
  # thresh : threshold(s) pi for active set selection
  # lambda.prop : threshold parameter theta for choosing a subset of the lambdas
  
  data.imp<-impute.fnc(data.obs,m)$data.imp
  
  # removing missing data occurring due to other reasons
  indices<-which(apply(data.imp,1,function(x) sum(is.na(x)))>0)
  data.imp<-data.imp[setdiff(1:nrow(data.imp),indices),]
  
  # model matrix
  datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                          paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp)
  
  # reordering columns from output to maintain the structure of the formula
  datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
  
  # running randomized LASSO to initiate a list of lambdas
  rand<-rbinom(2*(m+K)+1,1,0.5)
  rand.w<- 1/ifelse(rand==1,1,alph)
  p.fac<-rand.w
  wt<-rep(1,nrow(data.imp))
  if (!is.null(data.imp$weights)) wt<-1/data.imp$weights #inverse of sampling weights
  fit.randlas<-glmnet(as.matrix(datamat[,-1]),data.imp$Y,weights=wt,family='binomial',penalty.factor=p.fac)
  lambdas<-fit.randlas$lambda
  
  # running stable selection
  stable.sel<-array(0,c(2*m+2*K+2,length(lambdas),reps))
  for (i in 1:reps){
    print(i)
    # bootstrap sample
    data.obs.b<-data.obs[sample(1:nrow(data.obs),nrow(data.obs),replace=T),]
    
    # impute and run LASSO again
    data.imp<-impute.fnc(data.obs.b,m)$data.imp
    indices<-which(apply(data.imp,1,function(x) sum(is.na(x)))>0)
    data.imp<-data.imp[setdiff(1:nrow(data.imp),indices),]
    datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                            paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp)
    datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
    rand<-rbinom(2*(m+K)+1,1,0.5)
    rand.w<- 1/ifelse(rand==1,1,alph)
    p.fac<-rand.w
    wt<-rep(1,nrow(data.imp))
    if (!is.null(data.imp$weights)) wt<-1/data.imp$weights #inverse of sampling weights
    fit.randlas<-glmnet(as.matrix(datamat[,-1]),data.imp$Y,weights=wt,family='binomial',
                        penalty.factor=p.fac, lambda=lambdas)

    # matching lambdas with the output lambda list from glmnet; in some cases they can be different even though lambdas are initialized
    lambinds<-which(lambdas%in%fit.randlas$lambda)
    colinds<-sapply(1:length(fit.randlas$lambda),function(x) which(fit.randlas$lambda%in%lambdas[lambinds][x]))
    stable.sel[,lambinds,i]<-as.matrix(coef(fit.randlas)[,colinds])
    
  }
  
  ## removing lower values of lambda using lambda.prop
  stable.sel.ind<-(stable.sel!=0)*1
  prop.select<-apply(apply(stable.sel.ind,3,colSums),1,function(x) {k <- which(x!=0); mean(x[k])} )
  if (!is.null(lambda.prop)) {
    lam.valid.ind<-prop.select<=lambda.prop*nrow(stable.sel)
  } else lam.valid.ind<-1:length(lambdas)
  new.lambdas<-lambdas[lam.valid.ind]
  stable.sel.ind.omit<- stable.sel.ind[,lam.valid.ind,which(apply(stable.sel.ind,3,function(x) sum(x^2))!=0),drop=F]
  
  # some cleaning steps
  if (prod(dim(stable.sel.ind.omit))!=0) {
    
    lam.each.run<-apply(stable.sel.ind.omit,3, function(x) which(colSums(x)!=0))
    
    if (is.list(lam.each.run)) {
      stable.sel.ind.omit<-stable.sel.ind.omit[,,which(unlist(lapply(lam.each.run,length))!=0),drop=F]
    }
    
    if (prod(dim(stable.sel.ind.omit))!=0){  
      pi.j.lam<-apply(stable.sel.ind.omit,2,function(x){
        cols<-which(colSums(x)!=0)
        if (length(cols)==0){
          y<-rep(0,nrow(x))
        } else y<-rowMeans(x[,cols])
        return(y) 
      }) 
      pi.j.max<-apply(pi.j.lam,1,max)
      stable.set<-vector('list',length(thresh))
      for (t in 1:length(thresh)){
        stable.set[[t]]<-sort(which(pi.j.max>=thresh[t]))
        if(length(stable.set[[t]])==0) print(paste('no variables selected by Randomized LASSO for thresh =',thresh[t]))
      }
      
    } else {
      pi.j.max<-rep(0,2*(m+K+1))
      stable.set<-rep(0,0)
      stable.set.woint<-rep(0,0)
      prop.select<-NA
    }
  } else {
    pi.j.max<-rep(0,2*(m+K+1))
    stable.set<-rep(0,0)
    stable.set.woint<-rep(0,0)
    prop.select<-NA
  }
  
  
  ## now estimation steps
  coef.arr<-vector('list',length(thresh))
  for (t in 1:length(thresh)){
    coef.arr[[t]]<-matrix(0,reps2,(length(stable.set[[t]])))
  }
  
  for (j in 1:reps2){
    print(j)
    data.obs.b<-data.obs[sample(1:nrow(data.obs),nrow(data.obs),replace=T),]
    # impute and run LASSO again
    data.imp<-impute.fnc(data.obs.b,m)$data.imp
    indices<-which(apply(data.imp,1,function(x) sum(is.na(x)))>0)
    data.imp<-data.imp[setdiff(1:nrow(data.imp),indices),]
    wt<-rep(1,nrow(data.imp))
    if (!is.null(data.imp$weights)) wt<-1/data.imp$weights #inverse of sampling weights
    datamat<-model.matrix(as.formula(paste0('~Z+',paste0('S.',1:m,collapse='+'),'+',
                                            paste0('S.',1:m,'*Z',collapse='+'),'+W+W*Z')),data=data.imp)
    datamat<-datamat[,c(1,2,3:(m+2),(m+K+3):(2*m+K+2),(m+3):(m+K+2),(2*m+K+3):(2*m+2*K+2))]
    for (t in 1:length(thresh)){
      datamat.stable<-data.frame(datamat[,stable.set[[t]],drop=F])
      datamat.stable$Y<-data.imp$Y
      if(1%in%stable.set[[t]]) datamat.stable.g=datamat.stable[,-1,drop=F] else datamat.stable.g=datamat.stable
      fit.glm<-glm(Y~.,data=datamat.stable.g,weights=wt,family=binomial(link=link),control = list(maxit = 100))
      if (fit.glm$converged==T) coef.arr[[t]][j,]<-coef(fit.glm)
    }
  }
  beta.fin<-matrix(0,2*(m+K+1),length(thresh))
  beta.est<-vector('list',length(thresh))
  for (t in 1:length(thresh)){
    beta.est[[t]]<-colMeans(coef.arr[[t]][,colSums(coef.arr[[t]]^2)!=0 & !is.na(colMeans(coef.arr[[t]])),drop=F])
    if (length(beta.est[[t]])>0){
      beta.fin[(stable.set[[t]]),t]<-beta.est[[t]]
    } else beta.fin[,t]<-rep(NA,2*(m+K+1))
  }
  return(list(prop.select=prop.select,pi.j.max=pi.j.max,stable.set=stable.set,beta.las=beta.fin,beta.las.est=beta.est))
}


###############################################
## FUNCTIONS TO MEASURE MODEL PERFORMANCE #####
###############################################

# 7. To evaluate performance of a selected model in our simulation examples
predict.mod<-function(X,Y,beta,mod.pred='logit',mod.sim='logit',sig,mod.beta){
  if (mod.pred=='logit') {
    pred<-as.numeric(expit(data.matrix(cbind(1,X))%*%beta)>0.5)
    pred.prob<-expit(data.matrix(cbind(1,X))%*%beta)
  }
  else if (mod.pred=='probit') {
    pred<-as.numeric(pnorm(data.matrix(cbind(1,X))%*%beta)>0.5)
    pred.prob<-pnorm(data.matrix(cbind(1,X))%*%beta)
  }
  mod.err=NA
  if (!is.null(mod.beta)){
    if (mod.sim=='logit') {
      pred.mod<-expit(data.matrix(cbind(1,X))%*%mod.beta)
    }
    else if (mod.sim=='probit') {
      pred.mod<-pnorm(data.matrix(cbind(1,X))%*%mod.beta)
    }
    
    mod.err<-mean((pred.prob-pred.mod)^2)
  }
  misclass<-mean(pred!=Y)
  sens<-mean(pred[Y==1])
  spec<-1 - mean(pred[Y==0])
  s<-sum(beta[(m+3):(m+2+sig)]!=0)
  frac.s<-s/sig
  f<-sum(beta[(m+3+sig):(2*m+2)]!=0)
  frac.f<-f/(m-sig)
  return(list(pred=pred,pred.prob=pred.prob,pred.mod=pred.mod,scores=c(sens=sens,spec=spec,misclass=misclass,mod.err=mod.err,s=s,frac.s=frac.s,f=f,frac.f=frac.f)))
}


##############################################################
## BEST w*sens+(1-w)*spec with threshold limits
##############################################################

# 8. Optimal thresholds for best feature selection and/or prediction performance 
best.thresh<-function(mat,w.feat,w.perf,thresh.vec=thresh,thresh.min=0.5,thresh.max=0.9){
  thresh.n<-thresh.vec[which(thresh.vec>=thresh.min & thresh.vec<=thresh.max)]
  mat.n<-mat[,which(thresh.vec>=thresh.min & thresh.vec<=thresh.max)]
  score.vec1<-w.perf*mat.n[1,]+(1-w.perf)*mat.n[2,]
  score.vec2<-w.feat*mat.n[6,]+(1-w.feat)*(1-mat.n[8,])
  score.perf<-max(score.vec1,na.rm=T)
  score.feat<-max(score.vec2,na.rm=T)
  thresh.perf<-thresh.n[which(score.vec1==max(score.vec1,na.rm=T))]
  thresh.feat<-thresh.n[which(score.vec2==max(score.vec2,na.rm=T))]
  return(list(score.feat=score.feat,thresh.feat=thresh.feat,score.perf=score.perf,thresh.perf=thresh.perf))
}

# 9. feature selection and/or prediction performance for non BISS procedures
no.thresh<-function(mat,w.feat,w.perf){
  score.perf<-w.perf*mat[1]+(1-w.perf)*mat[2]
  score.feat<-w.feat*mat[6]+(1-w.feat)*(1-mat[8])
  return(list(score.feat=score.feat,score.perf=score.perf))
}

# 10. Optimal thresholds for least model error and/or misclassification error
best.misclassmod<-function(mat,thresh.vec=thresh,thresh.min=0.5,thresh.max=0.9){
  thresh.n<-thresh.vec[which(thresh.vec>=thresh.min & thresh.vec<=thresh.max)]
  mat.n<-mat[,which(thresh.vec>=thresh.min & thresh.vec<=thresh.max)]
  score.misclass<-min(mat.n[3,],na.rm=T)
  score.moderr<-min(mat.n[4,],na.rm=T)
  thresh.misclass<-thresh.n[which(mat.n[3,]==score.misclass)]
  thresh.moderr<-thresh.n[which(mat.n[4,]==score.moderr)]
  return(list(score.misclass=score.misclass,thresh.misclass=thresh.misclass,score.moderr=score.moderr,thresh.moderr=thresh.moderr))
}


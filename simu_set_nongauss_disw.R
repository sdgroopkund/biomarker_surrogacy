source(file="~/codes_submission/fun_general.R")

library(stats)
library(MASS)
library(mvtnorm)
library(mice)
library(glmnet)
require(e1071)
library(ncvreg)
library(rootSolve)

# 1. PARAMETERS I

p.IC<-1 
n<-500
n.t<-5000
m<-25
sig<-2
m1<-sig
sig2<-10
m2<-sig+sig2
K<-1 # dimension of W
script<-as.integer(commandArgs(trailingOnly=T)[1])

# 2. PARAMETERS II

val.w<-c(1,2,3)
prob.w<-rep(1/length(val.w),length(val.w))
lam.y1<-c(rep(1.5,m2),rep(1,m-m2))
lam.y0<-c(rep(1,m2),rep(1,m-m2))
alpha5<--0.6
alpha6<-0
beta00<-alpha5+sum(log(lam.y1/lam.y0))
beta01<-lam.y0-lam.y1
beta02<-alpha6
beta0<-c(beta00,beta01,beta02)
p<-rep(0,m+2)
p[1]<-0
p[2:(sig+1)]<-rep(10,sig)
p[m+2]<-0
beta10<-beta00-p[1]*abs(beta00)
beta11<-beta01-p[2:(m+1)]*abs(beta01)
beta12<-beta02-p[m+2]*abs(beta02)
beta1<-c(beta10,beta11,beta12)
beta<-c(beta00,beta10-beta00,beta01,beta11-beta01,beta02,beta12-beta02)

# 3. SIMULATING TRAINING DATA

U<-runif(n)
W_gen<-rep(0,n)
ind<-which(U<=prob.w[1])
W_gen[ind]<-val.w[1]
sum<-prob.w[1]
for (i in 2:(length(val.w)-1)){
  ind<-which(U> sum & (U<=sum+prob.w[i]))
  W_gen[ind]<-val.w[i]
  sum=sum+prob.w[i]
}
ind<-which(U>sum(prob.w[1:(length(val.w)-1)]))
W_gen[ind]<-val.w[length(val.w)]
Z_gen<-sample(c(rep(0,n/2),rep(1,n/2))) # treatment indicator
ind1<-which(Z_gen==0)
ind2<-which(Z_gen==1)
Y_star<-rbinom(n,1,expit(alpha5+alpha6*W_gen))
S_gen<- t(sapply(1:length(W_gen), function(t){
  if (Y_star[t]==1){
    rate.vec<-lam.y1*W_gen[t]
  }else{
    rate.vec<-lam.y0*W_gen[t]
  }
  rexp(m,rate.vec)
}))
X1<-cbind(rep(1,length(Y_star)),Z_gen,W_gen*S_gen,S_gen*(Z_gen*W_gen),W_gen,Z_gen*W_gen)
RR1<-X1%*%beta
Y1<-rbinom(length(Y_star),1,expit(RR1))
datain<-data.frame(S=S_gen,W=W_gen,Y=Y1,Z=Z_gen)
datain$IC<-((datain$Z==0) & (datain$Y==0))*rbinom(n,1,p.IC)

# 4. COLLECTING VARIABLES INTO DATA AND INDUCE MISSINGNESS

dataout<-datain
S<-as.matrix(subset(dataout,select=1:m))   ### S(1) for marker 1,..,m
Z<-dataout$Z
Y<-dataout$Y    ### disease outcome
W<-dataout$W    ### covariate
IC<-dataout$IC  ### indicator for placebo subjects who has S(1) measured
delta<-Z==1     ### indicator for validation subset
delta.CPV<-Z==1|IC==1  ### indicator for measurement of S(1)
dataout$delta<-delta
dataout$delta.CPV<-delta.CPV
data.obs<-within(dataout,{
  for (i in 1:m){
    eval(parse(text=paste0('S.',i,'[delta.CPV==F]<-NA')))
  }
})
data.obs<-data.obs[,-ncol(data.obs)]

# 5. SIMULATING TESTING DATA

U.t<-runif(n.t)
W_gen.t<-rep(0,n.t)
ind.t<-which(U.t<=prob.w[1])
W_gen.t[ind.t]<-val.w[1]
sum<-prob.w[1]
for (i in 2:(length(val.w)-1)){
  ind.t<-which(U.t> sum & (U.t<=sum+prob.w[i]))
  W_gen.t[ind.t]<-val.w[i]
  sum=sum+prob.w[i]
}
ind.t<-which(U.t>sum(prob.w[1:(length(val.w)-1)]))
W_gen.t[ind.t]<-val.w[length(val.w)]
Z_gen.t<-sample(c(rep(0,n.t/2),rep(1,n.t/2))) # treatment indicator
ind1.t<-which(Z_gen.t==0)
ind2.t<-which(Z_gen.t==1)
Y_gen.t<-rep(NA,n.t)
Y_star.t<-rbinom(n.t,1,expit(alpha5+alpha6*W_gen.t))
Y_z0.t<-Y_star.t[ind1.t]
Y_gen.t[ind1.t]<-Y_z0.t
S_gen.t<- t(sapply(1:length(W_gen.t), function(t){
  if (Y_star.t[t]==1){
    rate.vec.t<-lam.y1*W_gen.t[t]
  }else{
    rate.vec.t<-lam.y0*W_gen.t[t]
  }
  rexp(m,rate.vec.t)
}))
X1.t<-cbind(rep(1,length(Y_star.t)),Z_gen.t,W_gen.t*S_gen.t,S_gen.t*(Z_gen.t*W_gen.t),W_gen.t, Z_gen.t*W_gen.t)
RR1.t<-X1.t%*%beta
Y1.t<-rbinom(length(Y_star.t),1,expit(RR1.t))
datain.t<-data.frame(S=S_gen.t,W=W_gen.t,Y=Y1.t,Z=Z_gen.t)
datain.t$IC<-((datain.t$Z==0) & (datain.t$Y==0))*rbinom(n.t,1,p.IC)
datamat.t<-cbind(Z=datain.t$Z,SW=datain.t[,1:m]*datain.t$W,SZW=datain.t[,1:m]*datain.t$Z*datain.t$W,W=datain.t$W,WZ=datain.t$W*datain.t$Z)




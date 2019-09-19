source(file="~/codes_submission/fun_general.R")

library(stats)
library(MASS)
library(mvtnorm)
library(mice)
library(glmnet)
require(e1071)
library(ncvreg)

# 1. PARAMETERS I

p.IC<-1 
n<-500
n.t<-5000
m<-25
sig<-2
m1<-sig
sig2<-3
m2<-sig+sig2
K<-1 # dimension of W
script<-as.integer(commandArgs(trailingOnly=T)[1])

# 2. PARAMETERS II

mu.w<-0
sd.w<-1
alpha1<-c(rep(-1,m2),rep(0,m-m2))
alpha3<-c(rep(2,m2),rep(0,m-m2))
alpha2<-c(rep(0.5,m1),rep(0,m-m1)) 
alpha4<-alpha2
alpha5<--1.2
alpha6<-0
fac=0.5
sd.s.y1<-fac*rep(1,m)
sd.s.y0<-fac*rep(1,m)
rho.ss.y1<-0
rho.ss.y0<-0
rho.sw<-rep(0,m)
sigma.ss.y1<-diag(sd.s.y1)%*%((1-rho.ss.y1)*diag(m)+rho.ss.y1)%*%diag(sd.s.y1)
sigma.ss.y0<-diag(sd.s.y0)%*%((1-rho.ss.y0)*diag(m)+rho.ss.y0)%*%diag(sd.s.y0)
sigma<-sigma.ss.y1
beta00<-alpha5+0.5*(t(alpha3)%*%solve(sigma)%*%alpha3-t(alpha1)%*%solve(sigma)%*%alpha1)
beta01<-t(alpha1-alpha3)%*%solve(sigma)
beta02<-t(alpha3-alpha1)%*%solve(sigma)%*%alpha2+alpha6
beta0<-c(beta00,beta01,beta02)
p<-rep(0,m+2)
p[1]<-4
p[2:(sig+1)]<-rep(1,sig)
p[m+2]<-0
beta10<-beta00-p[1]*abs(beta00)
beta11<-beta01-p[2:(m+1)]*abs(beta01)
beta12<-beta02-p[m+2]*abs(beta02)
beta1<-c(beta10,beta11,beta12)
beta<-c(beta00,beta10-beta00,beta01,beta11-beta01,beta02,beta12-beta02)

# 3. SIMULATING TRAINING DATA

W_gen<-rnorm(n, mu.w,sd.w)
Z_gen<-sample(c(rep(0,n/2),rep(1,n/2))) # treatment indicator
ind<-which(Z_gen==0)
ind2<-which(Z_gen==1)
Y_star<-rbinom(n,1,expit(alpha5))
S_gen<- t(sapply(1:length(W_gen), function(t){
  if (Y_star[t]==1){
    mean<-alpha1+alpha2*W_gen[t]
    sigma<-sigma.ss.y1
  }else{
    mean<-alpha3+alpha2*W_gen[t]
    sigma<-sigma.ss.y0
  }
  rmvnorm(1,mean=mean,sigma=sigma)
}))
X1<-cbind(rep(1,length(Y_star)),Z_gen,S_gen,S_gen*Z_gen,W_gen, Z_gen*W_gen)
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

W_gen.t<-rnorm(n.t, mu.w,sd.w)
Z_gen.t<-sample(c(rep(0,n.t/2),rep(1,n.t/2))) # treatment indicator
ind.t<-which(Z_gen.t==0)
ind2.t<-which(Z_gen.t==1)
Y_gen.t<-rep(NA,n.t)
Y_star.t<-rbinom(n.t,1,expit(alpha5))
Y_z0.t<-Y_star.t[ind.t]
Y_gen.t[ind.t]<-Y_z0.t
S_gen.t<- t(sapply(1:length(W_gen.t), function(t){
  if (Y_star.t[t]==1){
    mean.t<-alpha1+alpha2*W_gen.t[t]
    sigma.t<-sigma.ss.y1
  }else{
    mean.t<-alpha3+alpha2*W_gen.t[t]
    sigma.t<-sigma.ss.y0
  }
  rmvnorm(1,mean=mean.t,sigma=sigma.t)
}))
X1.t<-cbind(rep(1,length(Y_star.t)),Z_gen.t,S_gen.t,S_gen.t*Z_gen.t,W_gen.t, Z_gen.t*W_gen.t)
RR1.t<-X1.t%*%beta
Y1.t<-rbinom(length(Y_star.t),1,expit(RR1.t))
datain.t<-data.frame(S=S_gen.t,W=W_gen.t,Y=Y1.t,Z=Z_gen.t)
datain.t$IC<-((datain.t$Z==0) & (datain.t$Y==0))*rbinom(n.t,1,p.IC)
datamat.t<-cbind(Z=datain.t$Z,S=datain.t[,1:m],SZ=datain.t[,1:m]*datain.t$Z,W=datain.t$W,WZ=datain.t$W*datain.t$Z)


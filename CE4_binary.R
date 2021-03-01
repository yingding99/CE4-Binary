library(multcomp)
library(msm)
library(geepack)
library(tableone)
library(dplyr)
library(Hmisc)
# Function for test
# Data in format: var1-"response"; var2-"Treatment"; var3-"Marker"
CE4_CI2 <- function(Data,alpha=0.05){
  # Data: contains 3 columns: Response (binary), Treatment ("Rx", "C"); Marker ("AA", "Aa", "aa")
  # calculate the population frequency of each marker category 
  Data<-Data[complete.cases(Data),]
  n0 <- nrow(Data[which(Data$Marker==0),])
  n1 <- nrow(Data[which(Data$Marker==1),])
  n2 <- nrow(Data[which(Data$Marker==2),])
  
  p.aa<-n2/(nrow(Data))
  p.Aa<-n1/(nrow(Data))
  p.AA<-1-p.aa-p.Aa
  
  # fit a loglinear regression model log(P(y=1)) = trt + marker + trt*marker
  fit.loglinear <- geeglm(Response ~ as.factor(Treatment) + as.factor(Marker) + as.factor(Treatment)*as.factor(Marker)+enrollage+smoke+SevScaleBL, data=Data,family=poisson(link="log"),id=SampleID,corstr  = "exchangeable")
  #summary(fit.loglinear)
  fit.LogCoeff<-fit.loglinear$coef
  fit.Logvar<-vcov(fit.loglinear)###variance-covariance matrix for beta_0 to beta_5
  mean.CE4<-rep(NA,4)
  mean.CE4[1]<-log(p.Aa*exp(fit.LogCoeff[3]+fit.LogCoeff[8])+p.aa*exp(fit.LogCoeff[4]+fit.LogCoeff[9]))-log(p.Aa*exp(fit.LogCoeff[3])+p.aa*exp(fit.LogCoeff[4]))
  mean.CE4[2]<-fit.LogCoeff[9]-log(p.AA+p.Aa*exp(fit.LogCoeff[3]+fit.LogCoeff[8]))+log(p.AA+p.Aa*exp(fit.LogCoeff[3]))
  mean.CE4[3]<-fit.LogCoeff[8]
  mean.CE4[4]<-fit.LogCoeff[9]-fit.LogCoeff[8]
  
  logr0<-fit.LogCoeff[2]
  logr1<-fit.LogCoeff[2]+fit.LogCoeff[8]
  logr2<-fit.LogCoeff[2]+fit.LogCoeff[9]
  logr12<-fit.LogCoeff[2]+log(p.Aa*exp(fit.LogCoeff[3]+fit.LogCoeff[8])+p.aa*exp(fit.LogCoeff[4]+fit.LogCoeff[9]))-log(p.Aa*exp(fit.LogCoeff[3])+p.aa*exp(fit.LogCoeff[4]))
  logr01<-fit.LogCoeff[2]+log(p.AA+p.Aa*exp(fit.LogCoeff[3]+fit.LogCoeff[8]))-log(p.AA+p.Aa*exp(fit.LogCoeff[3]))
  
  formr12 <- sprintf("~x2+(log(%s*exp(x3+x8)+ %s*exp(x4+x9))-log(%s*exp(x3)+ %s*exp(x4)))", p.Aa, p.aa, p.Aa, p.aa)
  formr01 <- sprintf("~(x2+log(%s+%s*exp(x3+x8))-log(%s+%s*exp(x3)))", p.AA, p.Aa, p.AA, p.Aa)
  logRR5.var <- deltamethod(list(~x2,~(x2+x8),~(x2+x9),as.formula(formr12), as.formula(formr01)), mean=fit.LogCoeff, cov=fit.Logvar, ses=FALSE)
  
  
  form1 <- sprintf("~(log(%s*exp(x3+x8)+ %s*exp(x4+x9))-log(%s*exp(x3)+ %s*exp(x4)))", p.Aa, p.aa, p.Aa, p.aa)
  form2 <- sprintf("~(x9-log(%s+%s*exp(x3+x8))+log(%s+%s*exp(x3)))", p.AA, p.Aa, p.AA, p.Aa)
  CE4.var <- deltamethod(list(as.formula(form1), as.formula(form2), ~x8, ~(x9-x8)), 
                         mean=fit.LogCoeff, cov=fit.Logvar, ses=FALSE)
  #CE4.var <- deltamethod(list(~(log(p.Aa*exp(x3+x5)+p.aa*exp(x4+x6))-log(p.Aa*exp(x3)+p.aa*exp(x4))),~(x6-log(p.AA+p.Aa*exp(x3+x5))+log(p.AA+p.Aa*exp(x3))),~x5,~(x6-x5)),mean=fit.LogCoeff,cov=fit.Logvar,ses=FALSE)
  
  q<-qmvnorm(1-alpha,tail=("both.tails"),corr=cov2cor(CE4.var))$quantile
  CI.CE.4<-matrix(rep(NA,8),nrow=4)
  CI.CE.4[,1]<-mean.CE4-q*sqrt(diag(CE4.var))
  CI.CE.4[,2]<-mean.CE4+q*sqrt(diag(CE4.var))
  colnames(CI.CE.4)<-c("lwr","upr")
  p.CE.4<-rep(NA,4)
  p.CE.4[1] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[1]/sqrt(diag(CE4.var)[1])),4),upper=rep(abs(mean.CE4[1]/sqrt(diag(CE4.var)[1])),4),corr=cov2cor(CE4.var))[1])
  p.CE.4[2] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[2]/sqrt(diag(CE4.var)[2])),4),upper=rep(abs(mean.CE4[2]/sqrt(diag(CE4.var)[2])),4),corr=cov2cor(CE4.var))[1])
  p.CE.4[3] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[3]/sqrt(diag(CE4.var)[3])),4),upper=rep(abs(mean.CE4[3]/sqrt(diag(CE4.var)[3])),4),corr=cov2cor(CE4.var))[1])
  p.CE.4[4] = 1-(pmvnorm(lower=rep(-abs(mean.CE4[4]/sqrt(diag(CE4.var)[4])),4),upper=rep(abs(mean.CE4[4]/sqrt(diag(CE4.var)[4])),4),corr=cov2cor(CE4.var))[1])
  p.CE4 = 1-(pmvnorm(lower=rep(-max(abs(mean.CE4/sqrt(diag(CE4.var)))),4),upper=rep(max(abs(mean.CE4/sqrt(diag(CE4.var)))),4),corr=cov2cor(CE4.var))[1])
  
  colnames(CI.CE.4)<-c("lwr","upr")
  
  return(list(p.CE.4=p.CE.4,CI.CE.4=CI.CE.4,est.CE.4=mean.CE4,p.CE4=p.CE4,logr0=logr0,logr1=logr1,logr2=logr2,logr01=logr01,logr12=logr12,logRR5.var=logRR5.var))
}





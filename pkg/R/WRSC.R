# these two C functions are from http://dornsife.usc.edu/assets/sites/239/docs/WRSC.txt


lintest_C<-function(x,y,regfun=tsreg_C,nboot=500,alpha=.05,xout=F,outfun=out,...){
#
# Test the hypothesis that the regression surface is a plane.
# Stute et al. (1998, JASA, 93, 141-149).
#
library(parallel)
set.seed(2)
if(identical(regfun,tshdreg))print('When using tshdreg, be sure to include RES=TRUE')
if(identical(regfun,Qreg))print('When using Qreg, be sure to include res.vals=TRUE')
if(identical(regfun,tsreg_C))print('When using C++, must use tshdregcpp with res.vals=TRUE')
if(identical(regfun,tshdreg_C))print('When using C++, must use tshdregcpp with res.vals=TRUE')
#  These last two commands create an error when using RStudio unless package WRScpp has been invoked.
x<-as.matrix(x)
d<-ncol(x)
temp<-elimna(cbind(x,y))
x<-temp[,1:d]
x<-as.matrix(x)
y<-temp[,d+1]
if(xout){
flag<-outfun(x)$keep
x<-x[flag,]
x<-as.matrix(x)
y<-y[flag]
}
mflag<-matrix(NA,nrow=length(y),ncol=length(y))
for (j in 1:length(y)){
for (k in 1:length(y)){
mflag[j,k]<-(sum(x[j,]<=x[k,])==ncol(x))
}
}
reg<-regfun(x,y,...)
yhat<-y-reg$residuals
print("Taking bootstrap sample, please wait.")
data<-matrix(runif(length(y)*nboot),nrow=nboot)
data<-sqrt(12)*(data-.5) # standardize the random numbers.
data=listm(t(data))
rvalb<-mclapply(data,lintests1,yhat,reg$residuals,mflag,x,regfun,mc.preschedule=TRUE,...)
# An n x nboot matrix of R values
rvalb=matl(rvalb)
rvalb<-rvalb/sqrt(length(y))
dstatb<-apply(abs(rvalb),2,max)
wstatb<-apply(rvalb^2,2,mean)
# compute test statistic
v<-c(rep(1,length(y)))
rval<-lintests1(v,yhat,reg$residuals,mflag,x,regfun,...)
rval<-rval/sqrt(length(y))
dstat<-max(abs(rval))
wstat<-mean(rval^2)
ib<-round(nboot*(1-alpha))
p.value.d<-1-sum(dstat>=dstatb)/nboot
p.value.w<-1-sum(wstat>=wstatb)/nboot
list(dstat=dstat,wstat=wstat,p.value.d=p.value.d,p.value.w=p.value.w)
}

scorsubMC<-function(isub,x,y,pr=FALSE,STAND=TRUE,corfun=corfun,cop=cop,CPP=TRUE,...){
isub=as.vector(isub)
if(!CPP)corbsub<-scor(x[isub],y[isub],plotit=FALSE,pr=FALSE,STAND=STAND,corfun=corfun,cop=cop,SEED=FALSE,...)$cor
if(CPP)corbsub<-scor_C(x[isub],y[isub],plotit=FALSE,pr=FALSE,STAND=STAND,corfun=corfun,cop=cop,SEED=FALSE,...)$cor
corbsub
}

scorci_C<-function(x,y,corfun=pcor,nboot=599,alpha=.05,SEED=TRUE,plotit=TRUE,
xlab = 'VAR 1', ylab ='VAR 2', STAND = TRUE, MGV = FALSE,...){
#
#   Compute a 1-alpha confidence interval for skipped correlation.
#   The default correlation is the Pearson's correlation.
#
#   The default number of bootstrap samples is nboot=599
#
library(WRScpp)
m1=cbind(x,y)
m1<-elimna(m1)  # Eliminate rows with missing values
nval=nrow(m1)
x<-m1[,1]
y<-m1[,2]
est<-scor_C(x,y,plotit=plotit,corfun=corfun,xlab=xlab,ylab=ylab,MGV=MGV,...)$cor      
if(SEED)set.seed(2) # set seed of random number generator so that                     
#             results can be duplicated.                                              
data<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)          
bvec<-apply(data,1,corbsub,x,y,scor_C,plotit=FALSE,...) # A 1 by nboot matrix.        
ihi<-floor((1-alpha/2)*nboot+.5)                                                      
ilow<-floor((alpha/2)*nboot+.5)                                                       
bsort<-sort(bvec)                                                                     
corci<-1                                                                              
corci[1]<-bsort[ilow]                                                                 
corci[2]<-bsort[ihi]                                                                  
phat <- sum(bvec < 0)/nboot                                                           
sig <- 2 * min(phat, 1 - phat)                                                        
list(cor.ci=corci,p.value=sig,cor.est=est)                                            
}                


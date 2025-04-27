
###############################################################################
#### The R code below conducts data fitting (for 1000 randomly selected subjects from the HP data) using the proposed methodology in the manuscript below. 
#### Manuscript: Fan Y and Bandyopadhyay D. (2025+). Regression-based individualized treatment rules for skewed bivariate outcomes with informative follow-up, (Submitted)
#### Authors: Yiwei Fan and Dipankar Bandyopadhyay
###############################################################################


rm(list = ls())

### Specify the directory you will be working from 

setwd("G:\\My Drive\\...\\")

source("myfunction.R")
mydf=read.csv("example_data.csv")
head(mydf)
dim(mydf)


################################ 
visiting model
#################################


install.packages("IrregLong", repos='http://cran.us.r-project.org')
library(IrregLong)

##### last-observation-carried-forward#######
#####covariates for the visiting model#########
mydf1=lagfn(mydf,lagvars=colnames(mydf)[c(3,4,6,9:14,18)],id="id",
            time="cumdelta",lagfirst=NA)
colnames(mydf1)
mydf1=mydf1[,c(16:28,5,7,8)]
mydf1=mydf1[,c(2,1,3,13,14:16,4:12)]
colnames(mydf1)
mydf1$CAL2.lag=(mydf1$CAL.lag)^2
mydf1$PD2.lag=(mydf1$PD.lag)^2

############################
##### Fiting the Cox model, stepwise feature selection and weight-fitting 
#############################

library(survival)
mph=coxph(Surv(cumdelta.lag,cumdelta,event)~.-id+cluster(id),
          data=mydf1)
coef(mph)
summary(mph)
####step-wise feature selection
step_mph=step(mph)
summary(step_mph)
output=data.frame(summary(step_mph)$coefficients)
output
####weight fitting
weight=iiw(mph,mydf1,"id","cumdelta",TRUE)
range(weight)


####################################
### Extract baseline CAL, PPD, and Treatment
####################################


first=match(unique(mydf$id),mydf$id)
mydf$BASE_CAL=rep(mydf[first,"CAL"],table(mydf$id))
mydf$BASE_PPD=rep(mydf[first,"PD"],table(mydf$id))
mydf$BASE_T=rep(mydf[first,"TX_IND1"],table(mydf$id))
colnames(mydf)



#########################################
### Fitting marginal models
#########################################


install.packages(c("glmtoolbox", "statmod", "fda", "tweedie"), repos='http://cran.us.r-project.org')
library(glmtoolbox)
library(statmod)
library(fda)
library(tweedie)

####CAL
mydf2=mydf[,c(3,5:14,17,19:21)]
last=c(first[-1]-1,nrow(mydf))
mydf2=mydf2[-last,]
colnames(mydf2)

#########time varying
library(fda)
range(mydf$cumdelta)
Bspline=create.bspline.basis(rangeval=c(0,95),
          norder=4,breaks=seq(0,95,length.out=4))
Bspline
temp=eval.basis(mydf$cumdelta[-first],Bspline)
n_B=length(Bspline$names)
subdf=mydf2
subdf=subdf[,c(12,1:5,10,13,15,11)]
subdf$CAL=mydf$CAL[-first]
colnames(subdf)

#######interation term
subdf$AT=subdf$AGE1*subdf$TX_IND1
subdf$BT=mydf2$BASE_CAL*subdf$TX_IND1
subdf$TT=mydf2$BASE_T*subdf$TX_IND1
subdf$ST1=mydf2$TOBACCO_IND1*subdf$TX_IND1
subdf$ST2=mydf2$TOBACCO_IND2*subdf$TX_IND1
subdf=as.matrix(subdf)
colnames(subdf)
r=2

######time-varing coefficients
i=4
tempX1=subdf[,i]*temp
i=10
tempX2=subdf[,i]*temp
subdf=cbind(subdf[,-c(4,10)],tempX1,tempX2)

#########################
### GEE for Tweedie distribution
##########################


library(glmtoolbox)
library(tweedie)
library(statmod)
Y1_wfit=glmgee(subdf[,2]~subdf[,-c(1,2)],
               id=subdf[,1],
               family=tweedie(var.power=r,link.power=0),
               weights=weight[-first])
print(QIC(Y1_wfit))
####estimated coeffcients
data.frame(summary(Y1_wfit)$coefficients)
Y1_wfit$coef

#########estimate phi
n_all=nrow(subdf)
a=sum((subdf[,2]-Y1_wfit$fitted.values)^2/(Y1_wfit$fitted.values)^r)
n_B=length(Bspline$names)
b=n_all-2-2*n_B
phihat1=a/b
phihat1


###############fitted values for \bar{Q}_1
newdf=data.frame(subdf[,-c(1,2)])
newdf$AT=0
newdf$BT=0
newdf$TT=0
newdf$ST1=0
newdf$ST2=0
newdf[,c(18:23)]=0
Y1_pred0=exp((c(Y1_wfit$coefficients)[-1])%*%t(newdf)+c(Y1_wfit$coefficients)[1])

newdf$AT=mydf2$AGE1
newdf$BT=mydf2$BASE_CAL
newdf$TT=mydf2$BASE_T
newdf$ST1=mydf2$TOBACCO_IND1
newdf$ST2=mydf2$TOBACCO_IND2
newdf[,c(18:23)]=temp
Y1_pred1=exp((c(Y1_wfit$coefficients)[-1])%*%t(newdf)+c(Y1_wfit$coefficients)[1])


##################################################
#####PPD
###############################################
colnames(mydf)
mydf3=mydf[-last,c(4:14,17,20,21)]

####time varying: AGE1, RACE1(TX_IND1)
subdf=mydf3[,-c(6,9)]
subdf=subdf[,c(10,1:5,8,11,12,9)]
subdf$PD=mydf$PD[-first]

####interation terms
subdf$AT=subdf$AGE1*subdf$TX_IND1
subdf$BT=subdf$BASE_PPD*subdf$TX_IND1
subdf$TT=mydf3$BASE_T*subdf$TX_IND1
subdf$ST1=mydf3$TOBACCO_IND1*subdf$TX_IND1
subdf$ST2=mydf3$TOBACCO_IND2*subdf$TX_IND1
subdf=as.matrix(subdf)
r=0

#####time-varying coefficients
i=4
tempX1=subdf[,i]*temp
i=10
tempX2=subdf[,i]*temp
subdf=cbind(subdf[,-c(4,10)],tempX1,tempX2)

####GEE for Tweedie distribution
Y2_wfit=glmgee(subdf[,2]~subdf[,-c(1,2)],
               id=subdf[,1],
               family=tweedie(var.power=r,link.power=0),
               weights=weight[-first])
print(QIC(Y2_wfit))
data.frame(summary(Y2_wfit)$coefficients)
Y2_wfit$coef

#########estimate phi
n_all=nrow(subdf)
a=sum((subdf[,2]-Y2_wfit$fitted.values)^2/(Y2_wfit$fitted.values)^r)
n_B=length(Bspline$names)
b=n_all-2-2*n_B
phihat2=a/b
phihat2


###############fitted values for \bar{Q}_2
newdf=data.frame(subdf[,-c(1,2)])
newdf$AT=0
newdf$BT=0
newdf$TT=0
newdf$ST1=0
newdf$ST2=0
newdf[,c(18:23)]=0
Y2_pred0=exp((c(Y2_wfit$coefficients)[-1])%*%t(newdf)+c(Y2_wfit$coefficients)[1])


newdf$AT=mydf3$AGE1
newdf$BT=mydf3$BASE_PPD
newdf$TT=mydf3$BASE_T
newdf$ST1=mydf3$TOBACCO_IND1
newdf$ST2=mydf3$TOBACCO_IND2
newdf[,c(18:23)]=temp
Y2_pred1=exp((c(Y2_wfit$coefficients)[-1])%*%t(newdf)+c(Y2_wfit$coefficients)[1])




############bivariate distribution regression################
tau1=0.4
tau2=0.4
####quantile-specific levels (covariate dependent but treatment independent)
tt1=qgamma(rep(tau1,n_all),shape=1/phihat1,scale=phihat1*Y1_pred0)
tt2=qgamma(rep(tau1,n_all),shape=1/phihat1,scale=phihat1*Y1_pred1)
tt3=qnorm(rep(tau2,n_all),mean=Y2_pred0,sd=sqrt(phihat2))
tt4=qnorm(rep(tau2,n_all),mean=Y2_pred1,sd=sqrt(phihat2))
tempy1=(tt1+tt2)/2
tempy2=(tt3+tt4)/2


######binary outcomes
OY1=(mydf2$CAL<=tempy1)
OY2=(mydf3$PD<=tempy2)
O_L=OY1*OY2
new1=pgamma(tempy1,shape=1/phihat1,scale=phihat1*Y1_wfit$fitted.values)
new2=pnorm(tempy2,mean=Y2_wfit$fitted.values,sd=sqrt(phihat2))


tempX=subdf[,-c(1,2,14:25)]
#########time-varying covariates
Bspline=create.bspline.basis(rangeval=c(0,95),
                             norder=4,breaks=seq(0,95,length.out=4))
temp=eval.basis(mydf$cumdelta[-first],Bspline)
n_B=length(Bspline$names)
tempX1=mydf3$AGE1*temp
tempX2=mydf3$TX_IND1*temp

#####covariates
tempX=cbind(tempX[,c(1:4)],mydf2$BASE_CAL,
            tempX[,c(5:7)],
            mydf2$BASE_CAL*mydf2$TX_IND1,
            tempX[,-c(1:7)],tempX1,tempX2)
tempX=cbind(1,tempX)
colnames(tempX)[c(6,10,11)]=c("BASE_CAL","BCT","BPT")
tempX=tempX[,c(1:2,9,13:ncol(tempX))]
colnames(tempX)


###########estimated results
set.seed(5)###
final=c()
est=optim(runif(ncol(tempX),-1,1),fn=ass_f,w=1/weight[-first],
          method="BFGS")
est$value
coef=round(est$par,3)
coef


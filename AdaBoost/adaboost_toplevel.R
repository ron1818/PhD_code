# adaboost toplevel
require(RSNNS)
require(MASS)
source('work/AdaBoost/myAdaBoostR2.R')
source('work/AdaBoost/myAdaBoostRT.R')
source('work/AdaBoost/myAdaBoostRPlus.R')
source('work/AdaBoost/myBEMBoost.R')
source('work/misc/pre_processing_fn.R')
# regression
pao<-read.table('work/RVFL/pao94.dat',header=T)

x<-my_scale_matrix(pao[,-5],c(0,1))$scale
y<-my_scale_matrix(pao[,5,drop=F],c(0,1))$scale
x.train<-x[1:30,]
y.train<-y[1:30]
x.test<-x[31:38,]
y.test<-y[31:38]

BEM.model<-BEMBoost(x.train,y.train,x.test,BEM=0.1)
R2.model<-AdaBoost.R2(x.train,y.train,x.test)
RT.model<-AdaBoost.RT(x.train,y.train,x.test)
RP.model<-AdaBoost.Plus(x.train,y.train,x.test)

plot(y.test, col='red',type='l')
lines(BEM.model$aggregated.predict,col='green')
lines(R2.model$aggregated.predict,col='yellow')
lines(as.numeric(RT.model$aggregated.predict),col='blue')
lines(RP.model$aggregated.predict,col='pink')

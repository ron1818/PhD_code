## load required package
library(forecast)
library(xts)
library(TSA)
library(e1071)
library(frbs)
library(fUnitRoots)

## set working directory

setwd("C:/Users/renye/Dropbox/MATLAB/wind")


## change zero speed to a very small amount
#speed.ts[speed.ts==0]<-1e-16

# # split by months
# speed.dec.ts=window(speed.ts, start=as.timeDate('2010-12-01 00:00'), end=as.timeDate('2010-12-31 23:50'))
# speed.jan.ts=window(speed.ts, start=as.timeDate('2011-01-01 00:00'), end=as.timeDate('2011-01-31 23:50'))
# speed.feb.ts=window(speed.ts, start=as.timeDate('2011-02-01 00:00'), end=as.timeDate('2011-02-28 23:50'))
# speed.mar.ts=window(speed.ts, start=as.timeDate('2011-03-01 00:00'), end=as.timeDate('2011-03-31 23:50'))
# speed.apr.ts=window(speed.ts, start=as.timeDate('2011-04-01 00:00'), end=as.timeDate('2011-04-30 23:50'))
# speed.may.ts=window(speed.ts, start=as.timeDate('2011-05-01 00:00'), end=as.timeDate('2011-05-31 23:50'))
# speed.jun.ts=window(speed.ts, start=as.timeDate('2011-06-01 00:00'), end=as.timeDate('2011-06-30 23:50'))
# speed.jul.ts=window(speed.ts, start=as.timeDate('2011-07-01 00:00'), end=as.timeDate('2011-07-31 23:50'))
# speed.aug.ts=window(speed.ts, start=as.timeDate('2011-08-01 00:00'), end=as.timeDate('2011-08-31 23:50'))
# speed.sep.ts=window(speed.ts, start=as.timeDate('2011-09-01 00:00'), end=as.timeDate('2011-09-30 23:50'))
# speed.oct.ts=window(speed.ts, start=as.timeDate('2011-10-01 00:00'), end=as.timeDate('2011-10-31 23:50'))
# speed.nov.ts=window(speed.ts, start=as.timeDate('2011-11-01 00:00'), end=as.timeDate('2011-11-30 23:50'))

# inspect TS
my_inspect<-function(x, month_name)
{
	# acf and pacf
	acf.test<-acf(x, main=paste(month_name, "ACF Plot"))
	pacf.test<-pacf(x, main=paste(month_name, "PACF Plot"))
	# eacf test
	
	# print('Extended ACF Test')
	eacf.test<-eacf(x)
	# unit root test
	# print('ADF Unit Root Test')
	urdf.test<-urdfTest(x)
	# print('ADF Unit Root Test, diff=1')
	urdff.test<-urdfTest(diff(x)[-1])
	return(list(acf.test<-acf.test, pacf.test<-pacf.test, eacf.test<-eacf.test, urdf.test=urdf.test, urdff.test=urdff.test))
}

# arima fit
my_arima_fit<-function(x,month_name,order)
{
	AIC.mat<-array(0, dim=c(4,4))
	rcount=1
	AR.range<-order[1]+(0:3)
	MA.range<-order[3]+(0:3)
	for (i in AR.range) 
	{ #AR
      ccount=1
		for (j in MA.range)
		{ #MA
			fit<-arima(x, order=c(i,order[2],j))
			AIC.mat[rcount,ccount]=fit$aic
			ccount=ccount+1
		}
		rcount=rcount+1
	}
	# show AIC table
	print(AIC.mat)
	# find min AIC
	min.AIC.index<-which(AIC.mat==min(AIC.mat), arr.ind=T)
	print(min.AIC.index)
	# fit best ARIMA
	fit<-arima(x, order=c(AR.range[min.AIC.index[1]],order[2],MA.range[min.AIC.index[2]]))
	# tsdiag(fit)
	# summary(fit)
	# ARCH Test
	# McLeod.Li.test(fit)
	return(list(fit=fit, aicmat=AIC.mat, aicidx=min.AIC.index))
}


# arima forecast
my_arima_forecast<-function(trn,tst,fit,HORIZON)
{
	combined=rbind(trn,tst)# combine trn and tst

	# delatst=max(length(fit$model$phi),HORIZON) # find matrn delatst
	pred<-array(0, dim=c(length(tst),5)) # initialize predicted version

	for (i in (1:length(tst))) # out of sample forecast
	{
		frac<-combined[i:(i+length(trn)-HORIZON)]
		fit2<-Arima(frac, model=fit)
		f<-forecast(fit2, h=HORIZON)
		# lower 80, 95, mean, upper 80, 95
		pred[i,]<-c(f$lower[HORIZON,], f$mean[HORIZON], f$upper[HORIZON,])
	}
	res<-tst-pred[,3] #mean
	MAE<-mean(abs(res))
	RMSE<-sqrt(mean(res^2))
	MAPE<-mean(abs(res/mean(tst)))*100
	# save testing_residuals, error measures, ARIMA parameters
	return(list(pred=pred, res=res, MAE=MAE, RMSE=RMSE, MAPE=MAPE))
}

# rearrange data to vector
my_vector<-function(trn,tst,HORIZON,F_COUNT)
{
	trn.length<-length(trn)-F_COUNT-HORIZON
	trn.data<-array(0,dim=c(trn.length, F_COUNT))
	trn.labels<-trn[(F_COUNT+HORIZON+1):length(trn)]
	for (i in (1:F_COUNT))
	{
		trn.data[,i]=trn[i:(i+trn.length-1)]
	}

        total=rbind(trn,tst)
        tst.data=array(0, dim=c(length(tst), F_COUNT))
        tst.labels=tst
        for (i in (1:F_COUNT))
	{
		tst.data[,i]=total[(i+length(trn)-F_COUNT-HORIZON)+(0:(length(tst)-1))]
	}

	return(list(trndata=trn.data, trnlabels=trn.labels, tstdata=tst.data, tstlabels=tst.labels))
}
	
# svr 
my_svr<-function(trn,tst,cost_range,gamma_range,epsilon_range)
{	
	MSE.mat<-array(0, dim=c(length(cost_range), length(gamma_range), length(epsilon_range)))
	count1=1
	# with 5-fold CV
	for (cost in cost_range)
	{
		count2=1		
		for (gamma in gamma_range)
		{
			count3=1			
			for (epsilon in epsilon_range)
			{
				m<-svm(trn$data, trn$labels, cost=cost, epsilon=epsilon, gamma=gamma, cross=5)
				# average MSE
				MSE.mat[count1,count2,count3]<-median(m$MSE)
				count3=count3+1
			}
			count2=count2+1
		}
		count1=count1+1
	}
	# select best performed m
	bestparam<-which(MSE.mat==min(MSE.mat),arr.ind=T)
	best.param<-c(cost_range[bestparam[1,1]], gamma_range[bestparam[1,2]], epsilon_range[bestparam[1,3]])
	# build best SVR
	m<-svm(trn$data, trn$labels, cost=best.param[1], gamma=best.param[2], epsilon=best.param[3])
	# test
	pred<-predict(m, tst$data)
	# error measures
	res<-tst$labels-pred #mean
	MAE<-mean(abs(res))
	RMSE<-sqrt(mean(res^2))
	MAPE<-mean(abs(res/mean(tst$labels)))*100
	return(list(model=m, MSE=MSE.mat, best.param=best.param, pred=pred, res=res, MAE=MAE, RMSE=RMSE, MAPE=MAPE))
}


# DENFIS
my_denfis<-function(trn,tst,cost_range,gamma_range,epsilon_range)
{	
	MSE.mat<-array(0, dim=c(length(cost_range), length(gamma_range), length(epsilon_range)))
	count1=1
	# with 5-fold CV
	for (cost in cost_range)
	{
		count2=1		
		for (gamma in gamma_range)
		{
			count3=1			
			for (epsilon in epsilon_range)
			{
				m<-svm(trn$data, trn$labels, cost=cost, epsilon=epsilon, gamma=gamma, cross=5)
				# average MSE
				MSE.mat[count1,count2,count3]<-median(m$MSE)
				count3=count3+1
			}
			count2=count2+1
		}
		count1=count1+1
	}
	# select best performed m
	bestparam<-which(MSE.mat==min(MSE.mat),arr.ind=T)
	best.param<-c(cost_range[bestparam[1,1]], gamma_range[bestparam[1,2]], epsilon_range[bestparam[1,3]])
	# build best SVR
	m<-svm(trn$data, trn$labels, cost=best.param[1], gamma=best.param[2], epsilon=best.param[3])
	# test
	pred<-predict(m, tst$data)
	# error measures
	res<-tst$labels-pred #mean
	MAE<-mean(abs(res))
	RMSE<-sqrt(mean(res^2))
	MAPE<-mean(abs(res/mean(tst$labels)))*100
	return(list(model=m, MSE=MSE.mat, best.param=best.param, pred=pred, res=res, MAE=MAE, RMSE=RMSE, MAPE=MAPE))
}

# main program
my_main<-function(ts, month, RATIO, HORIZON, F_COUNT)
{
	# Split into training and testing
	trn.ts<-ts[1:round(RATIO*length(ts))]
	tst.ts<-ts[(round(RATIO*length(ts))+1):length(ts)]
	# Fit ARIMA
	arima.fit<-auto.arima(as.numeric(trn.ts))
	# Forecast ARIMA
	arima.result<-my_arima_forecast(trn.ts, tst.ts, arima.fit, HORIZON)
	# Save ARIMA result
	# SVM regression
	svm.data<-my_vector(trn.ts, tst.ts, HORIZON, F_COUNT)
	svm.trn<-list(data=svm.data$trndata, labels=svm.data$trnlabels)
	svm.tst<-list(data=svm.data$tstdata, labels=svm.data$tstlabels)
	cost_range<-10**(-2:3)
	gamma_range<-(1/F_COUNT)*(1:F_COUNT)
	epsilon_range<-10**(-4:0)
	svm.result<-my_svr(svm.trn, svm.tst, cost_range, gamma_range, epsilon_range)

	# ARIMA-SVM regression
	arimasvm.data<-my_vector(arima.fit$residuals, arima.result$res, HORIZON[1], F_COUNT)
	arimasvm.trn<-list(data=arimasvm.data$trndata, labels=arimasvm.data$trnlabels)
	arimasvm.tst<-list(data=arimasvm.data$tstdata, labels=arimasvm.data$tstlabels)
	# Forecast
	arimasvm.nonlin<-my_svr(arimasvm.trn, arimasvm.tst, cost_range, gamma_range, epsilon_range)
	arimasvm.result<-arima.result$pred+arimasvm.nonlin$pred
	res<-tst.ts-arimasvm.result[,3]
	MAE<-mean(abs(res))
	RMSE<-sqrt(mean(res^2))
	MAPE<-mean(abs(res/mean(tst.ts)))*100
	arimasvm.result=list(fit=arimasvm.result, MAE=MAE, RMSE=RMSE, MAPE=MAPE) 
	
	# return(list(arima.coef=arima.fit$coef, arima.MAE=arima.result$MAE, arima.MAPE=arima.result$MAPE, arima.RMSE=arima.result$RMSE, svm.coef=svm.result$best.param, svm.MAE=svm.result$MAE, svm.MAPE=svm.result$MAPE, svm.RMSE=svm.result$RMSE, as.coef=arimasvm.nonlin$best.param, as.MAE=MAE, as.MAPE=MAPE, as.RMSE=RMSE))
	month.result<-list(month=month, HORIZON=HORIZON, arima.coef=arima.fit$coef, arima.MAE=arima.result$MAE, arima.MAPE=arima.result$MAPE, arima.RMSE=arima.result$RMSE, svm.coef=svm.result$best.param, svm.MAE=svm.result$MAE, svm.MAPE=svm.result$MAPE, svm.RMSE=svm.result$RMSE, as.coef=arimasvm.nonlin$best.param, as.MAE=MAE, as.MAPE=MAPE, as.RMSE=RMSE)
	# return(month.result)
	save(month.result, file=paste(month,'_', HORIZON, '_result.Rdata', sep=''))
}



## define month name
month=c('Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov')
## define trn:tst ratio
RATIO=0.7
## define HORIZON
HORIZON=c(6,12,24)
## define F_COUNT
F_COUNT=24

## open csv file 
wind<-read.csv("dataset/NDBC/41001h2011.csv")

## convert to ts format
speed.ts<-xts(wind[,3], as.timeDate(wind[,1]))
dir.ts<-xts(wind[,2], as.timeDate(wind[,1]))
speed.ts[speed.ts==0]<-0.1

## split (xts) by month
monthly.speed.ts<-split(speed.ts, f='months', drop=FALSE)
#notice it is a list variable, need [[]] to decompose
mapply(my_main, monthly.speed.ts, month, RATIO, HORIZON[1],F_COUNT)
mapply(my_main, monthly.speed.ts, month, RATIO, HORIZON[2],F_COUNT)
mapply(my_main, monthly.speed.ts, month, RATIO, HORIZON[3],F_COUNT)

# Tabulate result
SVM.coef.mat12=matrix(0, nrow=12, ncol=3) # c g e
ARIMASVM.coef.mat12=matrix(0, nrow=12, ncol=3) # c g e
error.mat12=matrix(0, nrow=36, ncol=3) # row:MAE, RMSE, MAPE, col: ARIMA, SVM, AS

for (i in 1:length(month))
{
	filename<-paste(month[i],'_', HORIZON[2], '_result.Rdata', sep='')
	load(filename)
	SVM.coef.mat12[i,]=month.result$svm.coef
	ARIMASVM.coef.mat12[i,]=month.result$as.coef
	error.mat12[3*i-2,]=c(month.result$arima.MAE, month.result$svm.MAE, month.result$as.MAE)#MAE
	error.mat12[3*i-1,]=c(month.result$arima.RMSE, month.result$svm.RMSE, month.result$as.RMSE)#MAE
	error.mat12[3*i,]=c(month.result$arima.MAPE, month.result$svm.MAPE, month.result$as.MAPE)#MAE
}

# Tabulate result
SVM.coef.mat24=matrix(0, nrow=12, ncol=3) # c g e
ARIMASVM.coef.mat24=matrix(0, nrow=12, ncol=3) # c g e
error.mat24=matrix(0, nrow=36, ncol=3) # row:MAE, RMSE, MAPE, col: ARIMA, SVM, AS

for (i in 1:length(month))
{
	filename<-paste(month[i],'_', HORIZON[3], '_result.Rdata', sep='')
	load(filename)
	SVM.coef.mat24[i,]=month.result$svm.coef
	ARIMASVM.coef.mat24[i,]=month.result$as.coef
	error.mat24[3*i-2,]=c(month.result$arima.MAE, month.result$svm.MAE, month.result$as.MAE)#MAE
	error.mat24[3*i-1,]=c(month.result$arima.RMSE, month.result$svm.RMSE, month.result$as.RMSE)#MAE
	error.mat24[3*i,]=c(month.result$arima.MAPE, month.result$svm.MAPE, month.result$as.MAPE)#MAE
}

ts2<-read.csv('2hour.csv')
MAE2=array(0,dim=c(3,3))
MAPE2=array(0,dim=c(3,3))
RMSE2=array(0,dim=c(3,3))
MAE2[1,3]<-wilcox.test(ts2$ARIMA.MAE, ts2$A.D.MAE, paired=T,alternative='g')$p.value
MAE2[2,3]<-wilcox.test(ts2$DENFIS.MAE, ts2$A.D.MAE, paired=T,alternative='g')$p.value
MAE2[3,2]<-wilcox.test(ts2$A.D.MAE, ts2$DENFIS.MAE, paired=T,alternative='g')$p.value
MAE2[1,2]<-wilcox.test(ts2$ARIMA.MAE, ts2$DENFIS.MAE, paired=T,alternative='g')$p.value
MAE2[2,1]<-wilcox.test(ts2$DENFIS.MAE, ts2$ARIMA.MAE, paired=T,alternative='g')$p.value
MAE2[3,1]<-wilcox.test(ts2$A.D.MAE, ts2$ARIMA.MAE, paired=T,alternative='g')$p.value
MAPE2[1,3]<-wilcox.test(ts2$ARIMA.MAPE, ts2$A.D.MAPE, paired=T,alternative='g')$p.value
MAPE2[2,3]<-wilcox.test(ts2$DENFIS.MAPE, ts2$A.D.MAPE, paired=T,alternative='g')$p.value
MAPE2[3,2]<-wilcox.test(ts2$A.D.MAPE, ts2$DENFIS.MAPE, paired=T,alternative='g')$p.value
MAPE2[1,2]<-wilcox.test(ts2$ARIMA.MAPE, ts2$DENFIS.MAPE, paired=T,alternative='g')$p.value
MAPE2[3,1]<-wilcox.test(ts2$A.D.MAPE, ts2$ARIMA.MAPE, paired=T,alternative='g')$p.value
MAPE2[2,1]<-wilcox.test(ts2$DENFIS.MAPE, ts2$ARIMA.MAPE, paired=T,alternative='g')$p.value
RMSE2[1,3]<-wilcox.test(ts2$ARIMA.RMSE, ts2$A.D.RMSE, paired=T,alternative='g')$p.value
RMSE2[2,3]<-wilcox.test(ts2$DENFIS.RMSE, ts2$A.D.RMSE, paired=T,alternative='g')$p.value
RMSE2[3,2]<-wilcox.test(ts2$ARIMA.RMSE, ts2$DENFIS.RMSE, paired=T,alternative='g')$p.value
RMSE2[1,2]<-wilcox.test(ts2$A.D.RMSE, ts2$DENFIS.RMSE, paired=T,alternative='g')$p.value
RMSE2[3,1]<-wilcox.test(ts2$A.D.RMSE, ts2$ARIMA.RMSE, paired=T,alternative='g')$p.value
RMSE2[2,1]<-wilcox.test(ts2$DENFIS.RMSE, ts2$ARIMA.RMSE, paired=T,alternative='g')$p.value

ts6<-read.csv('6hour.csv')
MAE6=array(0,dim=c(3,3))
MAPE6=array(0,dim=c(3,3))
RMSE6=array(0,dim=c(3,3))
MAE6[1,3]<-wilcox.test(ts6$ARIMA.MAE, ts6$A.D.MAE, paired=T,alternative='g')$p.value
MAE6[2,3]<-wilcox.test(ts6$DENFIS.MAE, ts6$A.D.MAE, paired=T,alternative='g')$p.value
MAE6[3,2]<-wilcox.test(ts6$A.D.MAE, ts6$DENFIS.MAE, paired=T,alternative='g')$p.value
MAE6[1,2]<-wilcox.test(ts6$ARIMA.MAE, ts6$DENFIS.MAE, paired=T,alternative='g')$p.value
MAE6[2,1]<-wilcox.test(ts6$DENFIS.MAE, ts6$ARIMA.MAE, paired=T,alternative='g')$p.value
MAE6[3,1]<-wilcox.test(ts6$A.D.MAE, ts6$ARIMA.MAE, paired=T,alternative='g')$p.value
MAPE6[1,3]<-wilcox.test(ts6$ARIMA.MAPE, ts6$A.D.MAPE, paired=T,alternative='g')$p.value
MAPE6[2,3]<-wilcox.test(ts6$DENFIS.MAPE, ts6$A.D.MAPE, paired=T,alternative='g')$p.value
MAPE6[3,2]<-wilcox.test(ts6$A.D.MAPE, ts6$DENFIS.MAPE, paired=T,alternative='g')$p.value
MAPE6[1,2]<-wilcox.test(ts6$ARIMA.MAPE, ts6$DENFIS.MAPE, paired=T,alternative='g')$p.value
MAPE6[3,1]<-wilcox.test(ts6$A.D.MAPE, ts6$ARIMA.MAPE, paired=T,alternative='g')$p.value
MAPE6[2,1]<-wilcox.test(ts6$DENFIS.MAPE, ts6$ARIMA.MAPE, paired=T,alternative='g')$p.value
RMSE6[1,3]<-wilcox.test(ts6$ARIMA.RMSE, ts6$A.D.RMSE, paired=T,alternative='g')$p.value
RMSE6[2,3]<-wilcox.test(ts6$DENFIS.RMSE, ts6$A.D.RMSE, paired=T,alternative='g')$p.value
RMSE6[3,2]<-wilcox.test(ts6$ARIMA.RMSE, ts6$DENFIS.RMSE, paired=T,alternative='g')$p.value
RMSE6[1,2]<-wilcox.test(ts6$A.D.RMSE, ts6$DENFIS.RMSE, paired=T,alternative='g')$p.value
RMSE6[3,1]<-wilcox.test(ts6$A.D.RMSE, ts6$ARIMA.RMSE, paired=T,alternative='g')$p.value
RMSE6[2,1]<-wilcox.test(ts6$DENFIS.RMSE, ts6$ARIMA.RMSE, paired=T,alternative='g')$p.value

ts12<-read.csv('12hour.csv')
MAE12=array(0,dim=c(3,3))
MAPE12=array(0,dim=c(3,3))
RMSE12=array(0,dim=c(3,3))
MAE12[1,3]<-wilcox.test(ts12$ARIMA.MAE, ts12$A.D.MAE, paired=T,alternative='g')$p.value
MAE12[2,3]<-wilcox.test(ts12$DENFIS.MAE, ts12$A.D.MAE, paired=T,alternative='g')$p.value
MAE12[3,2]<-wilcox.test(ts12$A.D.MAE, ts12$DENFIS.MAE, paired=T,alternative='g')$p.value
MAE12[1,2]<-wilcox.test(ts12$ARIMA.MAE, ts12$DENFIS.MAE, paired=T,alternative='g')$p.value
MAE12[2,1]<-wilcox.test(ts12$DENFIS.MAE, ts12$ARIMA.MAE, paired=T,alternative='g')$p.value
MAE12[3,1]<-wilcox.test(ts12$A.D.MAE, ts12$ARIMA.MAE, paired=T,alternative='g')$p.value
MAPE12[1,3]<-wilcox.test(ts12$ARIMA.MAPE, ts12$A.D.MAPE, paired=T,alternative='g')$p.value
MAPE12[2,3]<-wilcox.test(ts12$DENFIS.MAPE, ts12$A.D.MAPE, paired=T,alternative='g')$p.value
MAPE12[3,2]<-wilcox.test(ts12$A.D.MAPE, ts12$DENFIS.MAPE, paired=T,alternative='g')$p.value
MAPE12[1,2]<-wilcox.test(ts12$ARIMA.MAPE, ts12$DENFIS.MAPE, paired=T,alternative='g')$p.value
MAPE12[3,1]<-wilcox.test(ts12$A.D.MAPE, ts12$ARIMA.MAPE, paired=T,alternative='g')$p.value
MAPE12[2,1]<-wilcox.test(ts12$DENFIS.MAPE, ts12$ARIMA.MAPE, paired=T,alternative='g')$p.value
RMSE12[1,3]<-wilcox.test(ts12$ARIMA.RMSE, ts12$A.D.RMSE, paired=T,alternative='g')$p.value
RMSE12[2,3]<-wilcox.test(ts12$DENFIS.RMSE, ts12$A.D.RMSE, paired=T,alternative='g')$p.value
RMSE12[3,2]<-wilcox.test(ts12$ARIMA.RMSE, ts12$DENFIS.RMSE, paired=T,alternative='g')$p.value
RMSE12[1,2]<-wilcox.test(ts12$A.D.RMSE, ts12$DENFIS.RMSE, paired=T,alternative='g')$p.value
RMSE12[3,1]<-wilcox.test(ts12$A.D.RMSE, ts12$ARIMA.RMSE, paired=T,alternative='g')$p.value
RMSE12[2,1]<-wilcox.test(ts12$DENFIS.RMSE, ts12$ARIMA.RMSE, paired=T,alternative='g')$p.value

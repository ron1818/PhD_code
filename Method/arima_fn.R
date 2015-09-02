### inspect TS
#require(xts)
#require(forecast)
#require(fUnitRoots)
#require(TSA)
my_acf_pacf<-function(x)
{
  # acf and pacf
  acf.test<-acf(x, main=("ACF Plot"))
  pacf.test<-pacf(x, main=("PACF Plot"))
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

# function that do multi-step, re-estimation/non-re-estimation ARIMA modelling and forecasting
my_arima<-function(trn,tst,h=1,is.online=TRUE,model=NULL){
	# check data type
	if(is.list(trn)){# trn, tst are lists
		train.data<-c(trn$data[1,-ncol(trn$data)],trn$data[,ncol(trn$data)])
		test.data<-tst$data[,ncol(tst$data)]
	}else{ # not list, must be ts, xts or vector
		train.data=trn
		test.data=tst
	}
  # combine
  total<- c(train.data, test.data)#convert back to 1-d ts
  total.ts<-ts(total,frequency=1) #convert to ts
	n=length(test.data) # for online forecasting
	
  if(is.null(model)){# no pre-determined model, need to use auto arima
    fit <- auto.arima(train.data, parallel=T, num.cores=8,stepwise=F) # train with arima
    order <- arimaorder(fit) # save arima order for later
  }else{ # use previously generated model
    fit <- model
  }
	# save predict 
  fc.mean <- matrix(0, nrow=n, ncol=h) # initialize
  fc.lower <- matrix(0, nrow=n, ncol=h) # initialize
  fc.upper <- matrix(0, nrow=n, ncol=h) # initialize
  counter=1
  refit<-list()
  for(i in 1:n) # loop on test data
  {  
  	# roll to get new training data, append testing data one by one
    x <- window(total.ts, start=i+1, end=length(train.data) + i)
    # check if it is online
    if(!is.online){#not online, do not update model
        refit[[i]]<-fit # just copy
    }else if(is.numeric(is.online)){ # refit in blocks
      if(counter%%round(n/is.online)==0){ # reach block length
        refit[[i]]<-auto.arima(x, parallel=T, num.cores=8,stepwise=F)
      }else{ # did not reach block length, keep previous model
      	# 1st data, copy model, else copy previous model
        if(i==1){
          refit[[i]]<-fit
        }else{
          refit[[i]] <- refit[[i-1]]
        }
      }
      counter=counter+1
    }else{ # is.online==TRUE, stepwise refit
      refit[[i]] <- auto.arima(x, parallel=T, num.cores=8,stepwise=F)
    }
    # forecast for next h step
    pred<-forecast(Arima(x, model=refit[[i]]),h=h)
    fc.mean[i,] <- pred$mean # point
    fc.lower[i,] <- pred$lower[,'95%']# point
    fc.upper[i,] <- pred$upper[,'95%'] # point 
  }
  return(list(model=refit,predict=fc.mean,predict.lower=fc.lower,predict.upper=fc.upper))
}

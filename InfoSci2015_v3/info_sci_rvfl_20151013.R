source('info_sci_header_20151013.R')

# pre-allocation
rvfl.error<-rvfl.pacf.error<-rvfl.sar.error<-array(NA,dim=c(2,HORIZON_length,2,2,2,TRIAL,length(month)))
rvfl.best.n_h<-rvfl.pacf.best.n_h<-rvfl.sar.best.n_h<-array(NA,dim=c(2,2,2,TRIAL,length(month)))
dimnames(rvfl.best.n_h)<-dimnames(rvfl.pacf.best.n_h)<-dimnames(rvfl.sar.best.n_h)<-
    list(c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'),1:TRIAL,1:length(month))
dimnames(rvfl.error)<-dimnames(rvfl.pacf.error)<-dimnames(rvfl.sar.error)<-
    list(c('RMSE','MAE'),1:HORIZON_length,c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'),1:TRIAL,1:length(month))

rev.idx.pacf<-rev.idx.sar<-list()
for(mth in 1:12){
    print(paste('month=', mth))
    window.name<-paste(year,month[mth],sep='-')
    # scale
    tmp.load<-my_scale_ts(load.xts[window.name,'DMD'],c(0,1)) # scale to 0,1 for sigmoid
    load.data.ts<-tmp.load$scaled

    # determine periodicity
    N=96*2
    # determine arima order, p, q, P, Q, etc.
    m<-auto.arima(load.data.ts[1:N])
    arima.order<-arimaorder(m)
    seasonal.data<-arima.order[7] # need to check if it is NA
    pacf.data<-pacf(load.data.ts[1:N],lag=48*2,plot=FALSE)
    idx.pacf<-which(abs(pacf.data$acf)>1.96/sqrt(N))
    idx.sar<-c(1:arima.order[1],1:arima.order[4]+arima.order[7])
    if(is.na(seasonal.data)){# no seasonal
        LAG=max(idx.pacf)# LAG is max p
    }else{# seasonal, LAG is two period
        LAG=seasonal.data*2
    }
    rev.idx.pacf[[mth]]<-max(idx.pacf)-idx.pacf+1 # reverse idx for trn and tst data, for selective feature
    rev.idx.sar[[mth]]<-max(idx.sar)-idx.sar+1 # reverse idx for trn and tst data, for selective feature


    # split to 50%-50%, 24 hour ahead, need 48 hour data
    tmp<-my_vector_split(load.data.ts,HORIZON=max(HORIZON_range),LAG=LAG,RATIO=RATIO)

    for (trial in 1:TRIAL){
        print(paste('trial=',trial))
        seed=seed.array[trial]

        trn<-tmp$trn
        tst<-tmp$tst
        trn$labels<-trn$labels[,HORIZON_range]
        tst$labels<-tst$labels[,HORIZON_range]
        # RVFL to test
        print('full')

        n_h_range=seq(1,ncol(trn$data))
        model<-my.general.cv.rvfl(trn$data,trn$labels,n_h_range=n_h_range,cv=5,seed=seed,
                                  x.test=tst$data,y.test=tst$labels,
                                  input.bias=c(TRUE,FALSE),
                                  hidden.bias=c(TRUE,FALSE),
                                  direct.link=c(TRUE,FALSE),
                                  active.function='sigmoid')
        rvfl.error[,,,,,trial,mth]<-model$error[c(5,1),,,,]
        rvfl.best.n_h[,,,trial,mth]<-model$best.n_h


        # pacf feature selection
        trn.pacf<-list(data=trn$data[,rev.idx.pacf[[mth]]],labels=trn$labels)
        tst.pacf<-list(data=tst$data[,rev.idx.pacf[[mth]]],labels=tst$labels)

        print('pacf')
        # RVFL to test
        n_h=seq(1,96/2)
        rvfl.pacf.model<-my.general.cv.rvfl(trn.pacf$data,trn.pacf$labels,n_h_range=n_h_range,cv=5,seed=seed,
                                            x.test=tst.pacf$data,y.test=tst.pacf$labels,
                                            input.bias=c(TRUE,FALSE),
                                            hidden.bias=c(TRUE,FALSE),
                                            direct.link=c(TRUE,FALSE),
                                            is.scale=FALSE,
                                            active.function='sigmoid')
        rvfl.pacf.error[,,,,,trial,mth]<-rvfl.pacf.model$error[c(5,1),,,,]
        rvfl.pacf.best.n_h[,,,trial,mth]<-rvfl.pacf.model$best.n_h

        # sar feature selection
        trn.sar<-list(data=trn$data[,rev.idx.sar[[mth]]],labels=trn$labels)
        tst.sar<-list(data=tst$data[,rev.idx.sar[[mth]]],labels=tst$labels)

        print('sar')
        # RVFL to test
        n_h=seq(1,96/4)
        rvfl.sar.model<-my.general.cv.rvfl(trn.sar$data,trn.sar$labels,n_h_range=n_h,cv=5,seed=seed,
                                           x.test=tst.sar$data,y.test=tst.sar$labels,
                                           input.bias=c(TRUE,FALSE),
                                           hidden.bias=c(TRUE,FALSE),
                                           direct.link=c(TRUE,FALSE),
                                           active.function='sigmoid')
        rvfl.sar.error[,,,,,trial,mth]<-rvfl.sar.model$error[c(5,1),,,,]
        rvfl.sar.best.n_h[,,,trial,mth]<-rvfl.sar.model$best.n_h
    }
}

save(rvfl.error,
     rvfl.best.n_h,
     rvfl.pacf.error,
     rvfl.pacf.best.n_h,
     rvfl.sar.error,
     rvfl.sar.best.n_h,
     rev.idx.pacf,
     rev.idx.sar,
     file=paste('rvfl_result_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))





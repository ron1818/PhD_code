source('info_sci_header_20151013.R')

# pre-allocation
arima.model<-rf.model<-rf.pacf.model<-rf.sar.model<-ann.model<-ann.pacf.model<-ann.sar.model<-svr.model<-svr.pacf.model<-svr.sar.model<-list()

ann.error<-svr.error<-rf.error<-arima.error<-persistent.error<-array(NA,c(14,length(month),HORIZON_length))
ann.pacf.error<-svr.pacf.error<-rf.pacf.error<-array(NA,c(14,length(month),HORIZON_length))
ann.sar.error<-svr.sar.error<-rf.sar.error<-array(NA,c(14,length(month),HORIZON_length))

arima.param<-arima.order<-ann.n_h<-svr.param<-list()
ann.pacf.n_h<-svr.pacf.param<-list()
ann.sar.n_h<-svr.sar.param<-list()
for(mth in 1:12){
    ann.n_h[[mth]]<-svr.param[[mth]]<-list()
    ann.pacf.n_h[[mth]]<-svr.pacf.param[[mth]]<-list()
    ann.sar.n_h[[mth]]<-svr.sar.param[[mth]]<-list()
    print(paste('month=',mth))
    window.name<-paste(year,month[mth],sep='-')
    # scale
    tmp.load<-my_scale_ts(load.xts[window.name,'DMD'],c(0,1)) # scale to 0,1 for sigmoid
    load.data.ts<-tmp.load$scaled

    # determine periodicity
    N=96*2
    # determine arima order, p, q, P, Q, etc.
    m<-auto.arima(load.data.ts[1:N])
    arima.order[[mth]]<-arimaorder(m)
    arima.param[[mth]]<-m$coef
    seasonal.data<-arima.order[[mth]][7] # need to check if it is NA
    pacf.data<-pacf(load.data.ts[1:N],lag=48*2,plot=FALSE)
    idx.pacf<-which(abs(pacf.data$acf)>1.96/sqrt(N))
    idx.sar<-c(1:arima.order[[mth]][1],1:arima.order[[mth]][4]+arima.order[[mth]][7])
    if(is.na(seasonal.data)){# no seasonal
        LAG=max(idx.pacf)# LAG is max p
    }else{# seasonal, LAG is two period
        LAG=seasonal.data*2
    }
    rev.idx.pacf<-max(idx.pacf)-idx.pacf+1 # reverse idx for trn and tst data, for selective feature
    rev.idx.sar<-max(idx.sar)-idx.sar+1 # reverse idx for trn and tst data, for selective feature


    # split to 50%-50%, 24 hour ahead, need 48 hour data
    tmp<-my_vector_split(load.data.ts,HORIZON=max(HORIZON_range),LAG=LAG,RATIO=RATIO)

    trn<-tmp$trn
    tst<-tmp$tst

    # persistent
    persistent.predict<-array(NA,c(nrow(tst$data),HORIZON_length))
    for (i in seq_along(HORIZON_range))
        persistent.predict[,i]<-tst$data[,LAG-48+HORIZON_range[i]]

    # arima model
    arima.model<-m
    # train with arima.model on trn and tst
    rev.trn.data<-t(apply(trn$data,1,rev))
    rev.tst.data<-t(apply(tst$data,1,rev))
    arima.fit<-array(NA,c(nrow(rev.trn.data),HORIZON))
    arima.predict<-array(NA,c(nrow(rev.tst.data),HORIZON))
    for (j in 1:nrow(rev.trn.data))
        arima.fit[j,]<-forecast(Arima(rev.trn.data[j,],model=arima.model),h=HORIZON)$mean

    for (j in 1:nrow(rev.tst.data))
        arima.predict[j,]<-forecast(Arima(rev.tst.data[j,],model=arima.model),h=HORIZON)$mean

    arima.fit<-arima.fit[,HORIZON_range]
    arima.predict<-arima.predict[,HORIZON_range]

    print(paste('full'))
    for (i in seq_along(HORIZON_range)){


        h<-HORIZON_range[i]
        print(paste('h=',h))

        trn.h<-trn
        tst.h<-tst
        trn.h$labels<-trn.h$labels[,h]
        tst.h$labels<-tst.h$labels[,h]

        # predict with machine learning
        ann.model<-my_cv_ann(trn.h,tst.h,method='mlp',n_h_range=seq(96,192,24),is.ts=TRUE,cv.criteria='MSE')
        svr.model<-my_cv_svm(trn.h,tst.h,is.ts=TRUE,c_range=10^(-2:2), e_range=10^(-2:0),cv.criteria='MSE')
        rf.model<-my_cv_RF(trn.h,tst.h,is.ts=TRUE)

        ann.n_h[[mth]][[i]]<-ann.model$n_h
        svr.param[[mth]][[i]]<-svr.model$best.param
        ann.error[,mth,i]<-my_forecasting_measure(ann.model$predict,tst.h$labels)
        svr.error[,mth,i]<-my_forecasting_measure(svr.model$predict,tst.h$labels)
        rf.error[,mth,i]<-my_forecasting_measure(rf.model$predict,tst.h$labels)
    }

    persistent.error[,mth,]<-my_forecasting_measure(persistent.predict,tst$labels[,HORIZON_range])
    arima.error[,mth,]<-my_forecasting_measure(arima.predict,tst$labels[,HORIZON_range])


    # pacf feature selection
    trn.pacf<-list(data=trn$data[,rev.idx.pacf],labels=trn$labels)
    tst.pacf<-list(data=tst$data[,rev.idx.pacf],labels=tst$labels)

    print(paste('pacf'))
    for (i in seq_along(HORIZON_range)){

        h<-HORIZON_range[i]
        print(paste('h=',h))


        trn.h<-trn
        tst.h<-tst
        trn.h$labels<-trn.h$labels[,h]
        tst.h$labels<-tst.h$labels[,h]

        # predict with machine learning
        ann.model<-my_cv_ann(trn.h,tst.h,method='mlp',n_h_range=seq(96,192,24),is.ts=TRUE,cv.criteria='MSE')
        svr.model<-my_cv_svm(trn.h,tst.h,is.ts=TRUE,c_range=10^(-2:2), e_range=10^(-2:0),cv.criteria='MSE')
        rf.model<-my_cv_RF(trn.h,tst.h,is.ts=TRUE)

        ann.pacf.n_h[[mth]][[i]]<-ann.model$n_h
        svr.pacf.param[[mth]][[i]]<-svr.model$best.param
        ann.pacf.error[,mth,i]<-my_forecasting_measure(ann.model$predict,tst.h$labels)
        svr.pacf.error[,mth,i]<-my_forecasting_measure(svr.model$predict,tst.h$labels)
        rf.pacf.error[,mth,i]<-my_forecasting_measure(ann.model$predict,tst.h$labels)
    }


    # sar feature selection
    trn.sar<-list(data=trn$data[,rev.idx.sar],labels=trn$labels)
    tst.sar<-list(data=tst$data[,rev.idx.sar],labels=tst$labels)

    print(paste('sar'))
    for (i in seq_along(HORIZON_range)){

        h<-HORIZON_range[i]
        print(paste('h=',h))


        trn.h<-trn
        tst.h<-tst
        trn.h$labels<-trn.h$labels[,h]
        tst.h$labels<-tst.h$labels[,h]

        # predict with machine learning
        ann.model<-my_cv_ann(trn.h,tst.h,method='mlp',n_h_range=seq(96,192,24),is.ts=TRUE,cv.criteria='MSE')
        svr.model<-my_cv_svm(trn.h,tst.h,is.ts=TRUE,c_range=10^(-2:2), e_range=10^(-2:0),cv.criteria='MSE')
        rf.model<-my_cv_RF(trn.h,tst.h,is.ts=TRUE)

        ann.sar.n_h[[mth]][[i]]<-ann.model$n_h
        svr.sar.param[[mth]][[i]]<-svr.model$best.param
        ann.sar.error[,mth,i]<-my_forecasting_measure(ann.model$predict,tst.h$labels)
        svr.sar.error[,mth,i]<-my_forecasting_measure(svr.model$predict,tst.h$labels)
        rf.sar.error[,mth,i]<-my_forecasting_measure(ann.model$predict,tst.h$labels)
    }
}

save(arima.order,arima.param,
     ann.error,rf.error,svr.error,ann.n_h,svr.param,
     ann.pacf.error,rf.pacf.error,svr.pacf.error,ann.pacf.n_h,svr.pacf.param,
     ann.sar.error,rf.sar.error,svr.sar.error,ann.sar.n_h,svr.sar.param,
     file=paste('rvfl_result_ann_svr_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))

save(arima.error,persistent.error,file=paste('rvfl_result_arima_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))





# pre-allocation
rvfl.error<-array(NA,dim=c(2,HORIZON_length,2,2,2,TRIAL,length(month)))
rvfl.best.n_h<-array(NA,dim=c(2,2,2,TRIAL,length(month)))
dimnames(rvfl.best.n_h)<-
    list(c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'),1:TRIAL,1:length(month))
dimnames(rvfl.error)<-
    list(c('RMSE','MAE'),1:HORIZON_length,c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'),1:TRIAL,1:length(month))

rev.idx.pacf<-rev.idx.sar<-list()


arima.model<-rf.model<-ann.model<-svr.model<-list()

ann.error<-svr.error<-rf.error<-arima.error<-persistent.error<-array(NA,c(14,length(month),HORIZON_length))

arima.param<-arima.order<-ann.n_h<-svr.param<-list()

arima.time<-array(NA,c(12,2))
ann.time<-svr.time<-rf.time<-array(NA,c(12,2))
rvfl.time<-array(NA,c(12,2))

for(mth in 5){
    ann.n_h[[mth]]<-svr.param[[mth]]<-list()

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
    #   pacf.data<-pacf(load.data.ts[1:N],lag=48*2,plot=FALSE)
    #   idx.pacf<-which(abs(pacf.data$acf)>1.96/sqrt(N))
    #   idx.sar<-c(1:arima.order[1],1:arima.order[4]+arima.order[7])
    #   if(is.na(seasonal.data)){# no seasonal
    #     LAG=max(idx.pacf)# LAG is max p
    #   }else{# seasonal, LAG is two period
    #     LAG=seasonal.data*2
    #   }
    #   rev.idx.pacf[[mth]]<-max(idx.pacf)-idx.pacf+1 # reverse idx for trn and tst data, for selective feature
    #   rev.idx.sar[[mth]]<-max(idx.sar)-idx.sar+1 # reverse idx for trn and tst data, for selective feature
    #
    #
    # split to 50%-50%, 24 hour ahead, need 48 hour data
    tmp<-my_vector_split(load.data.ts,HORIZON=max(HORIZON_range),LAG=LAG,RATIO=RATIO)

    trn<-tmp$trn
    tst<-tmp$tst

    # persistent
    persistent.predict<-array(NA,c(nrow(tst$data),HORIZON_length))
    for (i in seq_along(HORIZON_range))
        persistent.predict[,i]<-tst$data[,LAG-48+HORIZON_range[i]]

    # arima model
    arima.trn.time.start<-proc.time()
    m1<-auto.arima(load.data.ts[1:N])
    arima.model<-m
    # train with arima.model on trn and tst
    rev.trn.data<-t(apply(trn$data,1,rev))
    rev.tst.data<-t(apply(tst$data,1,rev))
    arima.fit<-array(NA,c(nrow(rev.trn.data),HORIZON))
    arima.predict<-array(NA,c(nrow(rev.tst.data),HORIZON))

    for (j in 1:nrow(rev.trn.data))
        arima.fit[j,]<-forecast(Arima(rev.trn.data[j,],model=arima.model),h=HORIZON)$mean
    arima.trn.time<-proc.time()-arima.trn.time.start
    arima.tst.time.start<-proc.time()
    for (j in 1:nrow(rev.tst.data))
        arima.predict[j,]<-forecast(Arima(rev.tst.data[j,],model=arima.model),h=HORIZON)$mean
    arima.tst.time<-proc.time()-arima.tst.time.start

    arima.fit<-arima.fit[,HORIZON_range]
    arima.predict<-arima.predict[,HORIZON_range]

    arima.time[mth,]<-c(arima.trn.time[3],arima.tst.time[3])

    print(paste('full'))
    #for (i in seq_along(HORIZON_range)){
    i=1
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
    ann.time[mth,]<-c(ann.model$trn.time[3],ann.model$tst.time[3])
    svr.time[mth,]<-c(svr.model$trn.time[3],svr.model$tst.time[3])
    rf.time[mth,]<-c(rf.model$trn.time[3],rf.model$tst.time[3])
    #}

    persistent.error[,mth,]<-my_forecasting_measure(persistent.predict,tst$labels[,HORIZON_range])
    arima.error[,mth,]<-my_forecasting_measure(arima.predict,tst$labels[,HORIZON_range])

    #RVFL
    #for (trial in 1:TRIAL){
    trial=1
    print(paste('trial=',trial))
    seed=seed.array[trial]

    trn<-tmp$trn
    tst<-tmp$tst
    trn$labels<-trn$labels[,HORIZON_range,drop=F]
    tst$labels<-tst$labels[,HORIZON_range,drop=F]
    # RVFL to test
    print('full')

    n_h_range=seq(1,ncol(trn$data))
    model<-my.general.cv.rvfl(trn$data,trn$labels,n_h_range=n_h_range,cv=5,seed=seed,
                              x.test=tst$data,y.test=tst$labels,
                              input.bias=c(TRUE),
                              hidden.bias=c(FALSE),
                              direct.link=c(TRUE),
                              active.function='sigmoid')
    rvfl.error[,,,,,trial,mth]<-model$error[c(5,1),,,,]
    rvfl.best.n_h[,,,trial,mth]<-model$best.n_h
    rvfl.time[mth,]=c(model$trn.time[3],model$tst.time[3])
    #}
}

save(arima.time,rvfl.time,ann.time,svr.time,rf.time,
     file=paste('rvfl_result_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))


df.trn<-data.frame(arima=arima.time[,1],
                   ann=ann.time[,1],
                   svr=svr.time[,1],
                   rf=rf.time[,1],
                   rvfl=rvfl.time[,1]
)

df.tst<-data.frame(arima=arima.time[,2],
                   ann=ann.time[,2],
                   svr=svr.time[,2],
                   rf=rf.time[,2],
                   rvfl=rvfl.time[,2]
)

df.gg<-data.frame(
    method=rep(c('sARIMA','ANN','SVR','RF','RVFL'),each=12),
    trn=c(arima.time[,1],
          ann.time[,1],
          svr.time[,1],
          rf.time[,1],
          rvfl.time[,1]),
    tst=c(arima.time[,2],
          ann.time[,2],
          svr.time[,2],
          rf.time[,2],
          rvfl.time[,2]),
    dataset=factor(rep(1:12,5))
)

postscript('rvfl_trn_time_barplot.eps',width=8,height=5,family='Times',horizontal=FALSE)
ggplot(df.gg,aes(x=dataset,y=trn,fill=method))+geom_bar(stat="identity",position='dodge')+theme(legend.position="top")+labs(fill="") +xlab("Dataset")+ylab('Time (s)')
dev.off()


postscript('rvfl_tst_time_barplot.eps',width=8,height=5,family='Times',horizontal=FALSE)
ggplot(df.gg,aes(x=dataset,y=tst,fill=method))+geom_bar(stat="identity",position='dodge')+theme(legend.position="top")+labs(fill="") +xlab("Dataset")+ylab('Time (s)')
dev.off()

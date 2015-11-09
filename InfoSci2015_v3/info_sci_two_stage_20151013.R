# pre-allocation

FIR.error<-rvfl.1stage.error<-rvfl.2stage.error<-array(NA,dim=c(2,HORIZON_length,TRIAL,length(month)))

for(mth in 1:12){
    print(paste('month=', mth))
    window.name<-paste(year,month[mth],sep='-')
    # scale
    tmp.load<-my_scale_ts(load.xts[window.name,'DMD'],c(0,1)) # scale to 0,1 for sigmoid
    load.data.ts<-tmp.load$scaled
    LAG=96
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


        rvfl.model.1pass<-my.general.rvfl(trn$data,trn$labels,rvfl.best.n_h[1,2,1,trial,mth],x.test=tst$data,y.test=tst$labels,seed=seed,
                                          has.input.bias=FALSE,
                                          has.hidden.bias=FALSE,
                                          has.direct.link=TRUE,
                                          act.fn='sigmoid')
        rvfl.model.2pass<-my.general.rvfl(trn$data,trn$labels,rvfl.best.n_h[1,2,2,trial,mth],x.test=tst$data,y.test=tst$labels,seed=seed,
                                          has.input.bias=FALSE,
                                          has.hidden.bias=FALSE,
                                          has.direct.link=FALSE,
                                          act.fn='sigmoid')

        # FIR model based on least square
        rvfl.FIR.model<-FIR.model<-list(fit=array(NA,dim(trn$labels)),predict=array(NA,dim(tst$labels)))
        for (h in seq_along(HORIZON_range)){
            trn.df<-data.frame(trn$data,labels=trn$labels[,h])
            tst.df<-data.frame(tst$data,labels=tst$labels[,h])
            lm.model<-lm(labels~.,trn.df)
            FIR.model$predict[,h]<-lm.predict<-predict(lm.model,tst.df)
            FIR.model$fit[,h]<-lm.model$fitted.value

            trn2.df<-data.frame(cbind(lm.model$fitted.value,rvfl.model.2pass$fit[,h,1]),labels=trn$labels[,h])
            tst2.df<-data.frame(cbind(lm.predict,rvfl.model.2pass$predict[,h,1]),labels=tst$labels[,h])
            names(tst2.df)=names(trn2.df)
            lm2.model<-lm(labels~.,trn2.df)
            rvfl.FIR.model$fit[,h]<-lm2.model$fitted.value
            rvfl.FIR.model$predict[,h]<-predict(lm2.model,tst2.df)
        }

        # error measures
        rvfl.1stage.error[,,trial,mth]<-my_forecasting_measure(rvfl.model.1pass$predict[,,1],tst$labels)[c(5,1),]
        rvfl.2stage.error[,,trial,mth]<-my_forecasting_measure(rvfl.FIR.model$predict,tst$labels)[c(5,1),]
        FIR.error[,,trial,mth]<-my_forecasting_measure(FIR.model$predict,tst$labels)[c(5,1),]
    }
}


save(rvfl.1stage.error,rvfl.2stage.error,FIR.error,
     file=paste('rvfl12stage_result_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))



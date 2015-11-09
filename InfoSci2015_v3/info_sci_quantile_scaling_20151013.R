source('info_sci_header_20151013.R')

# pre-allocation
rvfl.error<-array(NA,dim=c(2,HORIZON_length,2,2,2,TRIAL,length(month)))
rvfl.best.n_h<-array(NA,dim=c(2,2,2,TRIAL,length(month)))
dimnames(rvfl.best.n_h)<-
    list(c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'),1:TRIAL,1:length(month))
dimnames(rvfl.error)<-
    list(c('RMSE','MAE'),1:HORIZON_length,c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'),1:TRIAL,1:length(month))
rvfl.time<-array(NA, c(12,2))

for(mth in 1:12){
    
    print(paste('month=', mth))
    window.name<-paste(year,month[mth],sep='-')
    # scale
    tmp.load<-my_scale_ts(load.xts[window.name,'DMD'],c(0,1)) # scale to 0,1 for sigmoid
    load.data.ts<-tmp.load$scaled
    
    LAG=96
    # split to 50%-50%, 24 hour ahead, need 48 hour data
    tmp<-my_vector_split(load.data.ts,HORIZON=max(HORIZON_range),LAG=LAG,RATIO=RATIO)
    
    trn<-tmp$trn
    tst<-tmp$tst
    
    #RVFL
    for (trial in 1:TRIAL){
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
                                  input.bias=c(TRUE,FALSE),
                                  hidden.bias=c(TRUE,FALSE),
                                  direct.link=c(TRUE,FALSE),
                                  is.scale=TRUE,scale.method='quantile',
                                  active.function='sigmoid')
        rvfl.error[,,,,,trial,mth]<-model$error[c(5,1),,,,]
        rvfl.best.n_h[,,,trial,mth]<-model$best.n_h
        rvfl.time[mth,]=c(model$trn.time[3],model$tst.time[3])
    }
}


save(rvfl.error,rvfl.best.n_h,rvfl.time,
     file=paste('rvfl_result_scope_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))

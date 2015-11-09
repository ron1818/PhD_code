source('info_sci_header_20151013.R')
TRIAL=1

# pre-allocation
mse.grid<-list()

for(mth in 1:12){
    mse.grid[[mth]]<-list()
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

        n_h_range=seq(1,ncol(trn$data)/4)
        mse.grid[[mth]][[trial]]<-my.general.cv.with.scale.fn(trn$data,trn$labels,n_h_range,seed=seed)

    }
}



# debug
for (i in 1:21){
    print(min(mse.grid[[1]][[1]][i,,]))
}





save(rvfl.error,rvfl.best.n_h,rvfl.time,
     file=paste('rvfl_result_scope_',format(Sys.time(),'%Y%m%d'),'.RData',sep=''))

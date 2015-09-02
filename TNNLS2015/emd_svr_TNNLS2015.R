# how to cite:
# @ARTICLE{Ren2015TNNLS,
#          author = {Ye Ren and P. N. Suganthan},
#          title = {A Novel Empirical Mode Decomposition With Support Vector Regression for Wind Speed Forecasting},
#          journal = {IEEE Trans. Neural Netw. Learn. Syst.},
#          year = {2015},
#          volume = {PP},
#          pages = {1},
#          number = {99}
# }


setwd('./TNNLS2015')
# load packages
require(forecast)
require(timeDate)
require(xts)
require(EMD)
require(RSNNS)
# load functions
source('../Method/mlp_fn.R')
source('../Method/svm_fn.R')
source('../EMD/')
source('../misc/pre_processing_fn.R')
source('../misc/post_processing_fn.R')
source('../misc/cross_validation_fn.R')

# 4) read dataset
wind<-read.csv('41004h2011.csv')
wind.xts<-xts(wind[,'WSPD'],timeDate(wind[,1]))
org.data.xts<-wind.xts['2011-01-01/2011-01-31']

# 5) define parameter
HORIZON=12
F_COUNT=LAG=24
RATIO=0.7
MAX.IMF=6
n_h_range=seq(LAG,2*LAG,6)

# 10) scale
data.xts<-my_scale_ts(org.data.xts,c(0,1))$scaled

# 20) split
tmp<-my_vector_split(data.xts,HORIZON=HORIZON,LAG=LAG,RATIO=RATIO)
trn<-tmp$trn
tst<-tmp$tst

# 30) EMD to extract data
emd.xts<-emd(as.numeric(data.xts),max.imf=MAX.IMF)
nimf<-emd.xts$nimf

# 20) partition into trn and tst
emd.trn<-emd.tst<-list()
emd.svr.predict<-emd.ann.predict<-array(NA,dim=c(nrow(tst$labels),ncol(tst$labels),nimf+1))
for (i in 1:(nimf+1)){ # for each IMF, +1 is for residue
	if(i>nimf){ # for residue
		tmp<-my_vector_split(emd.xts$residue,HORIZON=HORIZON,LAG=F_COUNT,RATIO=RATIO)
	}else{ # for imf
		tmp<-my_vector_split(emd.xts$imf[,i],HORIZON=HORIZON,LAG=F_COUNT,RATIO=RATIO)
	}
	emd.trn$data[[i]]<-tmp$trn$data
	emd.trn$labels[[i]]<-tmp$trn$labels
	emd.tst$data[[i]]<-tmp$tst$data
	emd.tst$labels[[i]]<-tmp$tst$labels
	
	# 30) cv ann to determine best k and d
	emd.ann.p.model<-my_cv_ann(tmp$trn,tmp$tst,n_h_range,is.ts=TRUE,cv.criteria='MSE')
	emd.svr.p.model<-my_cv_svm(tmp$trn,tmp$tst,is.ts=TRUE,cv.criteria='MSE')
	
	# 40) save predicted value
	emd.ann.predict[,,i]<-emd.ann.p.model$predict
	emd.svr.predict[,,i]<-emd.svr.p.model$predict
}
# 50) take sum of IMFs
emd.ann.p.predict<-array(NA,dim=c(nrow(tst$labels),HORIZON))
for (h in 1:HORIZON)
	emd.ann.p.predict[,h]<-apply(emd.ann.predict[,h,],1,sum)

# 60) calculate error
emd.ann.p.error<-my_forecasting_measure(emd.ann.p.predict,tst$labels)

# 51) take sum of IMFs
emd.svr.p.predict<-array(NA,dim=c(nrow(tst$labels),HORIZON))
for (h in 1:HORIZON)
	emd.svr.p.predict[,h]<-apply(emd.svr.predict[,h,],1,sum)

# 61) calculate error
emd.svr.p.error<-my_forecasting_measure(emd.svr.p.predict,tst$labels)

#########################################################
# 110) form a large feature set
emd.trn.m.data<-matrix(unlist(emd.trn$data),ncol=F_COUNT*(MAX.IMF+1),byrow=F)
emd.tst.m.data<-matrix(unlist(emd.tst$data),ncol=F_COUNT*(MAX.IMF+1),byrow=F)
# # 111) scale to (0,1)

# put into list
emd.trn.m<-list(data=emd.trn.m.data,labels=trn$labels)
emd.tst.m<-list(data=emd.tst.m.data,labels=tst$labels)

# 120) cv ann
emd.ann.m.model<-my_cv_ann(emd.trn.m,emd.tst.m,n_h_range,is.ts=TRUE)
emd.svr.m.model<-my_cv_svr(emd.trn.m,emd.tst.m,is.ts=TRUE)

# 130) save predicted value
emd.ann.m.predict<-emd.ann.m.model$predict
emd.svr.m.predict<-emd.svr.m.model$predict

# 140) calculate error
emd.ann.m.error<-my_forecasting_measure(emd.ann.m.predict,tst$labels)
emd.svr.m.error<-my_forecasting_measure(emd.svr.m.predict,tst$labels)

### how to cite ###
# @INPROCEEDINGS{Ren2014CIEL,
#                author = {Ye Ren and Xueheng Qiu and P. N. Suganthan},
#                title = {{EMD} based AdaBoost-{BPNN} Method for Wind Speed Forecasting},
#                booktitle = {Proc. IEEE Symposium on Computational Intelligence and Ensemble Learning
#                             (CIEL'14)},
#   year = {2014},
#   address = {Orlando, US},
#   month = dec
# }

### 0) initialize ###
# setwd('./CIEL2013') # current working 
# required packages, if not installed, please use:
# install.packages() function to install the packages
require(RSNNS)
require(MASS)
require(xts)
require(timeDate)
require(forecast)
require(EMD)
# source user defined functions
source('../AdaBoost/adaboost_fn.R')
source('../misc/pre_processing_fn.R')
source('../misc/post_processing_fn.R')
source('../misc/cross_validation_fn.R')
source('../Method/mlp_fn.R')

# 4) read dataset
wind<-read.csv('41004h2011.csv')
wind.xts<-xts(wind[,'WSPD'],timeDate(wind[,1]))
org.data.xts<-wind.xts['2011-01-01/2011-01-31']

# 5) define parameter
HORIZON=5
LAG=24
RATIO=0.7
MAX.IMF=5
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

# 40) partition into trn and tst
emd.trn<-emd.tst<-list()
imf.ada.predict<-imf.ann.predict<-array(NA,dim=c(nrow(tst$labels),ncol(tst$labels),nimf+1))
for (i in 1:(nimf+1)){ # for each IMF, +1 is for residue
	if(i>nimf){ # for residue
		tmp<-my_vector_split(emd.xts$residue,HORIZON=HORIZON,LAG=LAG,RATIO=RATIO)
	}else{ # for imf
		tmp<-my_vector_split(emd.xts$imf[,i],HORIZON=HORIZON,LAG=LAG,RATIO=RATIO)
	}
	emd.trn$data[[i]]<-tmp$trn$data
	emd.trn$labels[[i]]<-tmp$trn$labels
	emd.tst$data[[i]]<-tmp$tst$data
	emd.tst$labels[[i]]<-tmp$tst$labels
	
	# 50) cv knn to determine best k and d
	emd.ann.model<-my_cv_ann(tmp$trn, tmp$tst, n_h_range, method='mlp', is.ts=TRUE)
	imf.ann.predict[,,i]<-emd.ann.model$predict
	
	# 60) loop on HORIZON for adaboost training
	emd.ada.model<-list()
	for (h in 1:HORIZON){
		emd.ada.model[[h]]<-AdaBoost.R2(emd.trn$data[[i]],emd.trn$labels[[i]][,h],emd.tst$data[[i]],method='mlp',maxIter=10)
		imf.ada.predict[,h,i]<-emd.ada.model[[h]]$aggregated.predict
	}
	
}
# 70) take sum of IMFs
emd.ada.predict<-emd.ann.predict<-array(NA,dim=c(nrow(tst$labels),HORIZON))
for (h in 1:HORIZON){
	emd.ann.predict[,h]<-apply(imf.ann.predict[,h,],1,sum)
	emd.ada.predict[,h]<-apply(imf.ada.predict[,h,],1,sum)
}

# 80) error measure
emd.ann.error<-my_forecasting_measure(emd.ann.predict,tst$labels,tst$data[,LAG])
emd.ada.error<-my_forecasting_measure(emd.ada.predict,tst$labels,tst$data[,LAG])

# how to cite:
# @ARTICLE{Ren2014JPEE,
#          author = {Ye Ren and P. N. Suganthan},
#          title = {{EMD}-k{NN} models for wind speed forecasting},
#          journal = {Journal of Power and Energy Engineering},
#          year = {2014},
#          volume = {2},
#          pages = {176-185},
#          number = {4}
# }


setwd('./JPEE2014')
# load packages
require(forecast)
require(timeDate)
require(xts)
require(EMD)
# load functions
source('../method/knn_fn.R')
source('../misc/pre_processing_fn.R')
source('../misc/post_processing_fn.R')
source('../misc/cross_validation_fn.R')

# input data
# wind speed data from 78 Marine Dr. Singapore
# use WK 17 data as an example
marine.data<-read.csv('marine_drive_v4.csv')
marine.xts<-xts(marine.data[,'ult_mean_speed'],timeDate(marine.data[,1]))
data.xts<-marine.xts['2013-05-21/2013-05-27']

HORIZON_array=1:12
F_COUNT=24 # feature count
RATIO=0.7
MAX.IMF=6

# KNN param
k.range=c(3,9,15,21)
f.select=c(2,4,8,12,18,24)
f.deselect=F_COUNT-f.select
f.range<-array(TRUE,c(length(f.select),F_COUNT))
for (i in seq_along(f.select))
	f.range[i,]=c(rep(FALSE,f.deselect[i]),rep(TRUE,f.select[i]))

# 5) split data into trn and tst without emd
tmp<-my_vector_split(data.xts,HORIZON=max(HORIZON_array),LAG=F_COUNT,RATIO=RATIO)
trn<-tmp$trn
tst<-tmp$tst

# EMD-knn-P: parallel knn on each imf
# 10) EMD to extract data
emd.xts<-emd(as.numeric(data.xts),max.imf=MAX.IMF)
nimf<-emd.xts$nimf
# 20) partition into trn and tst
emd.trn<-emd.tst<-list()
emd.knn.predict<-array(NA,dim=c(nrow(tst$labels),ncol(tst$labels),nimf+1))
for (i in 1:(nimf+1)){ # for each IMF, +1 is for residue
	if(i>nimf){ # for residue
		tmp<-my_vector_split(emd.xts$residue,HORIZON=max(HORIZON_array),LAG=F_COUNT,RATIO=RATIO)
	}else{ # for imf
		tmp<-my_vector_split(emd.xts$imf[,i],HORIZON=max(HORIZON_array),LAG=F_COUNT,RATIO=RATIO)
	}
	emd.trn$data[[i]]<-tmp$trn$data
	emd.trn$labels[[i]]<-tmp$trn$labels
	emd.tst$data[[i]]<-tmp$tst$data
	emd.tst$labels[[i]]<-tmp$tst$labels
	
	# 30) cv knn to determine best k and d
	emd.knn.p.model<-my_knn(tmp$trn,tmp$tst,k_range=k.range,f_range=f.range,plot.rmse=F,is.ts=TRUE)
	
	# 40) save predicted value
	emd.knn.predict[,,i]<-emd.knn.p.model$predict
}
# 50) take sum of IMFs
emd.knn.p.predict<-array(NA,dim=c(nrow(tst$labels),max(HORIZON_array)))
for (h in 1:max(HORIZON_array)){
  emd.knn.p.predict[,h]<-apply(emd.knn.predict[,h,],1,sum)
}

# 60) calculate error
emd.knn.p.error<-my_forecasting_measure(emd.knn.p.predict,tst$labels)

#########################################################
# EMD-knn-M: multiple feature, feature extraction by emd
# 110) form a large feature set
emd.trn.m.data<-matrix(unlist(emd.trn$data),ncol=F_COUNT*(MAX.IMF+1),byrow=F)
emd.tst.m.data<-matrix(unlist(emd.tst$data),ncol=F_COUNT*(MAX.IMF+1),byrow=F)
# # 111) scale to (0,1)

# put into list
emd.trn.m<-list(data=emd.trn.m.data,labels=trn$labels)
emd.tst.m<-list(data=emd.tst.m.data,labels=tst$labels)

# 120) cv knn to determine best k and d
f.range.m<-matrix(rep(f.range,MAX.IMF+1),ncol=F_COUNT*(MAX.IMF+1),byrow=FALSE)
emd.knn.m.model<-my_knn(emd.trn.m,emd.tst.m,k_range=k.range,f_range=f.range.m,plot.rmse=F,is.ts=TRUE)

# 130) save predicted value
emd.knn.m.predict<-emd.knn.m.model$predict

# 140) calculate error
emd.knn.m.error<-my_forecasting_measure(emd.knn.m.predict,tst$labels)

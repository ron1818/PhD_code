BEMBoost<-function( x, y, x.test=NULL, method='mlp', maxIter=100, BEM=0.5,... ){
	# Feely 2000
	# Big Error Margin Boosting Algorithm
	# Given initialTrainingSet, initialDistribution and BEM
	# D0[i] = initialDistribution[i] for all i
	# t = 0; errCount = 0 
	# 1. Populate currentTrainingSet using Dt
	# 2. Construct network ht and train it using currentTrainingSet.
	# 3. For each member of initialTrainingSet:
	# 	Present that element to the trained network and record the prediction error.
	# If the resulting prediction error is greater then BEM,
	# mark that element as incorrectly predicted.
	# errCount = errCount + 1
	# Otherwise, mark that element as correctly predicted.
	# 4. Calculate the UpFactor and DownFactor as follows:
	# 	errCount
	# initialTrainingSet UpFactor =
	# 	UpFactor
	# DownFactor 1 =
	# 	5. Cycle through each member of the originalTrainingSet and update each element’s
	# distribution in the next training set as follows:
	# 	If the element was incorrectly predicted by the current network:
	# 	distribution [] [] i distribution i UpFactor t+1 = t ×
	# If the element was correctly predicted by the current network:
	# 	distribution [] [] i distribution i DownFactor t+1 = t ×
	# 6. Populate a new training set according to the distribution just calculated.
	# 7. Repeat steps 1-6 until a given maximum number of networks has been created. 
	
	m<-nrow(x)
	n<-ncol(x)
	y<-matrix(y,nrow=m)
	
	T=maxIter # max iter
	#initialize weight distr and error
	epsilon<-beta<-rep(NA,T)
	upfactor<-downfactor<-rep(NA,T)
	fit<-e<-D<-array(NA,dim=c(T,m))
	D[1,]=rep(1/m, m)
	model<-list()
	# iterate
	errorcnt<-rep(0,T)
	
	# iterate
	t=1
	while (t<T){
		# bootstrap
		idx<-sample(m,replace=TRUE,prob=D[t,])
		tth_x=x[idx,]
		tth_y=y[idx,]
		# call weak learner
		model[[t]]=eval(parse(text=method))(tth_x,tth_y,...)
		
		# getback to hypothesis
		fit[t,]<-y_fit<-predict(model[[t]],x)
		AE=abs(y_fit-y) # absolute error
		# calculate error count AE>BEM
		is_AE_gt_BEM=AE>BEM
		errorcnt[t]=sum(is_AE_gt_BEM)
		# check if errorcnt is zero
		if(errorcnt[t]==0)
			break
		# calculate upfactor and downfactor
		upfactor[t]=m/errorcnt[t]
		downfactor[t]=1/upfactor[t]
		
		# update distr
		is_AE_le_BEM=!is_AE_gt_BEM
		
		D[t+1,is_AE_gt_BEM]=D[t,is_AE_gt_BEM]*upfactor[t]
		D[t+1,is_AE_le_BEM]=D[t,is_AE_le_BEM]*downfactor[t]
		# normalize
		D[t+1,]=D[t+1,]/sum(D[t+1,])
		t=t+1
	}
	max_t=t-1
	# final model:
	
	if(max_t>1)
		aggregated.fit<-apply(fit[1:max_t,],2,mean)
	else
		aggregated.fit<-fit[1:max_t,]
	
	# predict
	if(!is.null(x.test)){ # have test dataset
		m.test<-nrow(x.test)
		predict<-array(NA,dim=c(max_t,m.test))
		for (t in 1:max_t){
			predict[t,]<-predict(model[[t]],x.test)
		}
		
		if(max_t>1)
			aggregated.predict<-apply(predict,2,mean)
		else
			aggregated.fit<-predict
		
	}else{
		aggregated.predict=NULL
	}
	return(list(aggregated.fit=aggregated.fit,aggregated.predict=aggregated.predict,
							base.fit=fit,base.predict=predict,max.t=max_t))
}



#AdaBoost function for regression

AdaBoost.R2<-function ( x, y, x.test=NULL, method='mlp', loss.function='linear', maxIter=100, ...){
	# Boosting for Regression Transfer
	# David Pardoe and Peter Stone
	m<-nrow(x)
	n<-ncol(x)
	y<-matrix(y,nrow=m)
	
	T=maxIter # max iter
	#initialize weight distr and error
	epsilon<-beta<-D<-rep(NA,T)
	fit<-e<-w<-array(NA,dim=c(T,m))
	w[1,]=rep(1/m, m)
	model<-list()
	# iterate
	t=1
	while (t<T){
		idx<-sample(m,replace=TRUE,prob=w[t,])
		tth_x=x[idx,]
		tth_y=y[idx,]
		
		# call weak learner, CART
		model[[t]]=eval(parse(text=method))(tth_x,tth_y,...)
		# getback to hypothesis
		fit[t,]<-y_fit<-predict(model[[t]],x)
		# adjusted error for each instance
		loss<-abs(y-y_fit)
		D[t]=max(loss)
		
		if(loss.function=='linear')
			e[t,]=loss/D[t]
		else if(loss.function=='square')
			e[t,]=(loss/D[t])^2
		else # exponential
			e[t,]=1-exp(-loss/D[t])
		
		# adjusted error of hypothesis
		epsilon[t]<-sum(e[t,]*w[t,])
		if(epsilon[t]>=0.5)
			break
		# calculate beta
		beta[t]<-epsilon[t]/(1-epsilon[t])
		# update weight vector
		w[t+1,]=w[t,]*beta[t]^(1-e[t,])
		# normalize
		w[t+1,]=w[t+1,]/sum(w[t+1,])
		t=t+1
	}
	max_t=t-1
	
	# final model: weighted median
	ww=log(1/beta[1:max_t])
	
	aggregated.fit<-rep(NA,m)
	for(i in 1:m){
		aggregated.fit[i]<-weighted.median(fit[1:max_t,i],ww)
	}
	
	# predict
	if(!is.null(x.test)){ # have test dataset
		m.test<-nrow(x.test)
		predict<-array(NA,dim=c(max_t,m.test))
		for (t in 1:max_t){
			predict[t,]<-predict(model[[t]],x.test)
		}
		aggregated.predict<-rep(NA,m.test)
		for(i in 1:m.test){
			aggregated.predict[i]<-weighted.median(predict[,i],ww)
		}
		
	}else{
		aggregated.predict=NULL
	}
	return(list(aggregated.fit=aggregated.fit,aggregated.predict=aggregated.predict,
							base.fit=fit,base.predict=predict,beta=beta[1:max_t],weight=ww,max.t=max_t))
}



weighted.median<-function(x,w,...){
	n<-length(x)
	min.w<-min(w)
	expand.factor<-ceiling(1/min.w)
	w.int<-ceiling(expand.factor*w)
	cumsum.w.int<-cumsum(c(0,w.int))
	
	tmp<-rep(NA,sum(w.int))
	for ( i in 1:n){
		tmp[(cumsum.w.int[i]+1):cumsum.w.int[i+1]]<-rep(x[i],w.int[i])
	}
	
	return(median(tmp,...))
	
}

AdaBoost.RT<-function ( x, y, x.test=NULL, method='mlp', phi=0.5, maxIter=100, ...){
	# AdaBoost.RT: a Boosting Algorithm for Regression Problems
	#	D. P. Solomatine, D. L. Shrestha
	
	m<-nrow(x)
	n<-ncol(x)
	y<-matrix(y,nrow=m)
	
	T=maxIter # max iter
	#initialize weight distr and error
	epsilon<-beta<-D<-rep(NA,T)
	fit<-e<-w<-array(NA,dim=c(T,m))
	w[1,]=rep(1/m, m)
	
	# iterate
	t=1
	while (t<T){
		idx<-sample(m,replace=TRUE,prob=w[t,])
		tth_x=x[idx,]
		tth_y=y[idx,]
		
		# call weak learner, CART
		model[[t]]=eval(parse(text=method))(tth_x,tth_y,...)
		# getback to hypothesis
		fit[t,]<-y_fit<-predict(model[[t]],x)
		# adjusted error for each instance
		ARE<-abs((y-y_fit)/y)
		is_ARE_gt_phi=ARE>phi
		is_ARE_le_phi=ARE<=phi
		
		# calculate error of hypothesis
		epsilon[t]=sum(w[t,is_ARE_gt_phi])
		
		# calculate beta
		beta[t]<-epsilon[t]^2
		# update weight vector
		w[t+1,is_ARE_le_phi]=w[t,is_ARE_le_phi]*beta[t]
		w[t+1,is_ARE_gt_phi]=w[t,is_ARE_gt_phi]
		# normalize
		w[t+1,]=w[t+1,]/sum(w[t+1,])
		t=t+1
	}
	max_t=t-1
	
	# final model
	ww=log(1/beta[1:max_t])
	
	aggregated.fit<-ww%*%fit[1:max_t,]/sum(ww)
	
	# predict
	if(!is.null(x.test)){ # have test dataset
		m.test<-nrow(x.test)
		predict<-array(NA,dim=c(max_t,m.test))
		for (t in 1:max_t){
			predict[t,]<-predict(model[[t]],x.test)
		}
		aggregated.predict<-ww%*%predict/sum(ww)
	}else{
		aggregated.predict=NULL
	}
	return(list(aggregated.fit=aggregated.fit,aggregated.predict=aggregated.predict,
							base.fit=fit,base.predict=predict,beta=beta[1:max_t],weight=ww,max.t=max_t))
}


AdaBoost.Plus<-function ( x, y, x.test=NULL, method='mlp', phi=0.5, sigma=0.1, maxIter=100, ...){
	# ADABOOST+: an ensemble learning approach for estimating
	# weather-related outages in distribution systems
	# Padmavathy Kankanala, Sanjoy Das, and Anil Pahwa 
	
	#require(MASS)
	
	m<-nrow(x)
	n<-ncol(x)
	y<-matrix(y,nrow=m)
	
	T=maxIter # max iter
	#initialize weight distr and error
	epsilon<-beta<-D<-rep(NA,T)
	fit<-e<-w<-array(NA,dim=c(T,m))
	w[1,]=rep(1/m, m)
	model<-list()
	
	# iterate
	t=1
	while (t<T){
		idx<-sample(m,replace=TRUE,prob=w[t,])
		tth_x=x[idx,]
		tth_y=y[idx,]
		
		# call weak learner, CART
		model[[t]]=eval(parse(text=method))(tth_x,tth_y,...)
		# getback to hypothesis
		fit[t,]<-y_fit<-predict(model[[t]],x)
		# adjusted error for each instance
		ARE<-abs((y-y_fit)/y)
		is_ARE_gt_phi=ARE>phi
		is_ARE_le_phi=ARE<=phi
		
		# calculate error of hypothesis
		epsilon[t]=sum(w[t,is_ARE_gt_phi])
		
		# calculate beta
		beta[t]<-epsilon[t]^2
		# update weight vector
		w[t+1,is_ARE_le_phi]=w[t,is_ARE_le_phi]*beta[t]
		w[t+1,is_ARE_gt_phi]=w[t,is_ARE_gt_phi]
		# normalize
		w[t+1,]=w[t+1,]/sum(w[t+1,])
		t=t+1
	}
	max_t=t-1
	
	
	# the above is the same as adaboost.RT
	
	Y<-fit[1:max_t,]
	if(is.null(sigma))# no regularization
		delta<-t(ginv(Y))%*%y
	else
		delta<-t(ginv(t(Y)%*%Y+sigma*diag(m))%*%t(Y))%*%y
	
	aggregated.fit<-t(Y)%*%delta
	
	# predict
	if(!is.null(x.test)){ # have test dataset
		m.test<-nrow(x.test)
		predict<-array(NA,dim=c(max_t,m.test))
		for (t in 1:max_t){
			predict[t,]<-predict(model[[t]],x.test)
		}
		aggregated.predict<-t(predict)%*%delta
	}else{
		aggregated.predict=NULL
	}
	return(list(aggregated.fit=aggregated.fit,aggregated.predict=aggregated.predict,
							base.fit=fit,base.predict=predict,beta=beta[1:max_t],delta=delta,max.t=max_t))
}


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


AdaBoost.R2<-function ( x, y, x.test=NULL, method='mlp', loss.function='linear', maxIter=100, ...){
	# Boosting for Regression Transfer
	# David Pardoe and Peter Stone
	m<-nrow(x)
	n<-ncol(x)
	y<-matrix(y,nrow=m)
	#n_y<-ncol(y) # multi output case
	
	T=maxIter # max iter
	#initialize weight distr and error
	epsilon<-beta<-D<-rep(NA,T)
	e<-w<-array(NA,dim=c(T,m))
	fit<-array(NA,dim=c(T,m))
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

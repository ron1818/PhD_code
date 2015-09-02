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
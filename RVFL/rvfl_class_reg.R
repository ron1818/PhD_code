# define function sigmoid
sigmoid<-function(x){
	x<-as.matrix(x)
	# output is from 0 to 1
	return(1/(1+exp(-x)))
}

my_RVFL<-function(trn,tst=NULL,n_h_range=c(10,20),
									type=if(is.factor(trn$labels)) 'classification' else 'regression',
									is.roll=FALSE,is.reest=FALSE,is.ts=FALSE,is.online=FALSE,
									cv=5,cv.criteria='MSE'){
	
	# scale to [0,1] to fit sigmo
	scaled.trn.data<-my_scale_matrix(trn$data,c(0,1))	
	scaled.tst.data<-my_scale_matrix(tst$data,c(0,1),scaled.trn.data$range)	
	if(type=='regression'){
	scaled.trn.labels<-my_scale_matrix(trn$labels,c(0,1))	
	scaled.tst.labels<-my_scale_matrix(tst$labels,c(0,1),scaled.trn.labels$range)	
	}else{
		scaled.trn.labels<-trn$labels
		scaled.tst.labels<-tst$labels
	}
	
	
	# check if it is classification or regression
	if(type=='classification'){
		cv.criteria='F_score'
		is.roll=FALSE
	}

  trn.time.start<-proc.time()
	if(cv>1){
		# parallel training to get cv accuracy, # for ts data
		cv.list<-my_cv_partition(list(data=scaled.trn.data$scaled,labels=trn$labels),is.ts=is.ts,cv=cv)
		measure.grid<-array(0,dim=c(length(n_h_range),cv)) # initialize
		for(i in seq_along(n_h_range)) # loop to choose hidden layer number
		{
			for (j in (1:cv.list$k)){ # loop for k fold
				trn.data<-cv.list$trn.data[[j]]
				trn.labels<-cv.list$trn.labels[[j]]
				val.data<-cv.list$val.data[[j]]
				val.labels<-cv.list$val.labels[[j]]
				if(!is.vector(trn.labels)&&!is.factor(trn.labels)){
					trn.labels=trn.labels[,1]
					val.labels=val.labels[,1]
				}
				m<-RVFL(trn.data, trn.labels, n_h=n_h_range[i], x.test=val.data, y.test=val.labels)
				pred<-m$predict
				if(type!='classification'){
					measure.grid[i,j]<-my_forecasting_measure(pred,val.labels)[cv.criteria,1]
				}else{
					# reciprocal, small is better
					measure.grid[i,j]<-1/my_classification_measure(pred,val.labels)[[cv.criteria]]
				}
			}
		}
		
		# average out to find cv accuracy
		measure.score<-rowMeans(measure.grid,na.rm=TRUE)
		
		if (all(is.na(measure.score))) { # all are nan
			return(message('grid search cannot find a proper parameter set, please add class weight'))
		}
		
		# find best param from the grid
		best.param.idx<-which(measure.score==min(measure.score,na.rm=T),arr.ind=T)
		best.n_h<-n_h_range[best.param.idx]
		
	}else{# no need cv
		measure.grid<-best.param.idx<-NULL
		best.n_h<-n_h_range[1]
	}
	trn.time<-proc.time()-trn.time.start

  tst.time.start<-proc.time()
	best.m<-RVFL(scaled.trn.data$scaled, trn$labels, n_h=best.n_h, x.test=tst$data, y.test=tst$labels)
	tst.time<-proc.time()-tst.time.start
	return(list(model=best.m,predict=best.m$predict,best.n_h=best.n_h, trn.time=trn.time,tst.time=tst.time))
	
}


RVFL<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL,seed=NULL){
	# single output RVFL, if needs multiple output, call it at higher level
	# check y
	if(is.factor(y)){	
		type='classification'
	}else{
		type='regression'
	}

	# dim information  
	N<-nrow(x) # sample size
	n_i<-ncol(x) #feature size
	n_o=length(y)
	
	# initialize train parameter  
	# initialize test parameter  
	if(!is.null(x.test)){
		y.pred<-array(NA,dim=c(nrow(x.test),n_o))
	}else{
		y.pred<-NULL
	}
	
	# create random weights
	if(!is.null(seed)){ # with pre defined seed
		set.seed(seed)
	}
	tmp<-runif((n_h+1)*(n_i+1), -1, 1) # +1 for bias
	A<-matrix(tmp,nrow=n_h+1,ncol=n_i+1)

	
			x.new<-cbind(x,rep(1,N)) # has input bias
			RV<-A[1:n_h,]

		enhancement.node<-sigmoid(RV%*%t(x.new))

			enhancement.node.new<-t(enhancement.node) # no hidden bias,transpose
			all.node<-cbind(x.new,enhancement.node.new)

		# output layer Least square
	if(type=='classification') {		
		W<-ginv(all.node)%*%(c(-1,1)[y]) # keast square
		y.fit<-sign(all.node%*%W)
		y.fit<-factor(y.fit,c(-1,1))
	}else{
		W<-ginv(all.node)%*%y # keast square
		y.fit<-all.node%*%W
	}
	
	### testing
	if(!is.null(x.test)){
		N.t<-nrow(x.test)
		# has input.bias
		x.t.new<-cbind(x.test,rep(1,N.t))
		RV.t<-A[1:n_h,]
		
		en.t<-sigmoid(RV.t%*%t(x.t.new))
		# no hidden bias
		en.t.new<-t(en.t)
		# direct link
		all.node.t<-cbind(x.t.new,en.t.new)
		
		if(type=='classification') {
			y.pred<-sign(all.node.t%*%W)
			y.pred<-factor(y.pred,c(-1,1))
		}else{
			y.pred<-all.node.t%*%W
		}
	}	
		
	return(list(W=W,A=A,fit=y.fit,predict=y.pred))
}


my_bagging_RVFL<-function(trn,tst=NULL,B=100,...){
	# one pass to get best.n_h
	one.pass.m<-my_RVFL(trn,tst,...)
	best.n_h<-one.pass.m$best.n_h
	
	# bootstrap
	require(bootstrap)
	idx.length<-nrow(trn$data)
	boot<-bootstrap(1:idx.length, B, as.numeric)
	idx.boot<-boot$thetastar # oob.idx=-idx.boot
	
	model.list<-list()
	length(model.list)<-B
	single.pred<-array(0, dim=c(nrow(tst$data),B))
	oob.error<-single.error<-list()
	n_h_array<-rep(best.n_h,B)

	# for each bag
	for (j in (1:B)) # change to parallel
	{
		trn.bag<-trn$data[idx.boot[,j],] # training data
		trn.oob<-trn$data[-idx.boot[,j],] # oob data
		trn.bag.labels<-trn$labels[idx.boot[,j]] # trining labels
		trn.oob.labels<-trn$labels[-idx.boot[,j]] # oob labels
		m<-RVFL(trn.bag,trn.bag.labels,n_h=n_h_array[j], x.test=trn.oob,y.test=trn.oob.labels )
		model.list[[j]]<-m
		oob.pred<-m$predict  
		
		oob.error[[j]]<-my_classification_measure(oob.pred, trn.oob.labels)
		
		A<-m$A
		W<-m$W
		N.t<-nrow(tst$data)
		# has input.bias
		x.t.new<-cbind(tst$data,rep(1,N.t))
		RV.t<-A[1:n_h_array[j],]
		
		en.t<-sigmoid(RV.t%*%t(x.t.new))
		# no hidden bias
		en.t.new<-t(en.t)
		# direct link
		all.node.t<-cbind(x.t.new,en.t.new)
		
		y.pred<-sign(all.node.t%*%W)
		single.pred[,j]<-y.pred

		single.error[[j]]<-my_classification_measure(as.factor(single.pred[,j]), tst$labels)
	}
	bagging.pred<-my_majority_vote(single.pred)
	
	bagging.error<-my_classification_measure(bagging.pred$value,tst$labels)
	
# 	single.f_score=rep(NA,B)
# 	for(i in 1:B)
# 		single.f_score[i]<-single.error[[i]]$F_score
# 	summary(single.f_score)	
# 
# 	for(i in 1:B)
# 		print(oob.error[[i]]$F_score)
	
	return(list(single.predict=single.pred,bagging.predict=bagging.pred,single.error=single.error,bagging.error=bagging.error))
}

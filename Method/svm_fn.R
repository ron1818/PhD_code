### svr with e1071 
require(e1071)

my_svm<-function(x,y,x.test=NULL,...){
  model<-svm(x=x,y=y,...)
  if(!is.null(x.test))
    y.predict<-predict(model,x.test)
  else
    y.predict<-NULL
  
  return(list(model=model,predict=y.predict))
}

# parallel processing k-fold CV to train SVR
my_cv_svm<-function(trn,tst=NULL,kernel='radial',
								 type=if(is.factor(trn$labels)) 'C-classification' else 'eps-regression',
								 c_range=10^(-4:4),
								 g_range=if (is.vector(trn$data)) 1 else 1/ncol(trn$data),
								 e_range=10^(-4:0),
								 d_range=2:5,
								 is.roll=FALSE,is.reest=FALSE,
								 is.ts=FALSE,is.online=FALSE,
								 cv=5,cv.criteria='MSE',...
) {
	# is.reest is to retrain svr with new horizon or just keep on using the original model
	# is.online is to retrain svr with new data or just keep old
	# is.roll is to roll forecast or direct forecast
	# x <- list(data,labels)
	# or x<- list(trn.data, trn.labels,val.data,val.labels, k*)
	
	#   require(e1071)
	#   require(bootstrap)
	#   require(foreach)
	if(type=='C-classification'){
		e_range=0.1
		cv.criteria='F_score'
		is.roll=FALSE
	}
	
	# SVR parameter flatten to an array, ease to parallel
	if (kernel == 'linear')#linear
	{#linear
		# preallocation for parallel
		score_grid_length=length(c_range)*length(e_range)
		param_c=rep(0, score_grid_length) # parameter, c
		param_e=rep(0, score_grid_length) # parameter, e
		param_d=rep(0, score_grid_length) # dummy
		param_g=rep(0, score_grid_length) # dummy
		
		coarse_score_grid=rep(0, score_grid_length) # accuracy, c e
		counter=1
		for (c in 1:length(c_range)){
			for (e in 1:length(e_range)){
				param_c[counter]=c_range[c]
				param_e[counter]=e_range[e]
				counter=counter+1
			}
		}
		
		
	}else if (kernel == 'polynomial')
	{#poly
		# preallocation for parallel
		score_grid_length=length(c_range)*length(e_range)*length(d_range)*length(e_range)
		param_c=rep(0, score_grid_length) # parameter, c
		param_e=rep(0, score_grid_length) # parameter, e
		param_d=rep(0, score_grid_length) # parameter, d
		param_g=rep(0, score_grid_length) # parameter, g
		
		coarse_score_grid=rep(0, score_grid_length) # accuracy, c e
		counter=1
		for (c in 1:length(c_range)){
			for (g in 1:length(g_range)){
				for (e in 1:length(e_range)){
					for (d in 1:length(d_range)){
						param_c[counter]=c_range[c]
						param_g[counter]=g_range[g]
						param_e[counter]=e_range[e]
						param_d[counter]=d_range[d]
						counter=counter+1
					}
				}
			}
		}
	}else
	{ # rbf or sigm
		# preallocation for parallel
		score_grid_length=length(c_range)*length(g_range)*length(e_range)
		param_c=rep(0, score_grid_length) # parameter, c
		param_e=rep(0, score_grid_length) # parameter, e
		param_d=rep(0, score_grid_length) # dummy
		param_g=rep(0, score_grid_length) # parameter, g
		
		coarse_score_grid=rep(0, score_grid_length) # accuracy, c e
		counter=1
		for (c in 1:length(c_range)){
			for (g in 1:length(g_range)){
				for (e in 1:length(e_range)){
					param_c[counter]=c_range[c]
					param_g[counter]=g_range[g]
					param_e[counter]=e_range[e]
					counter=counter+1
				}
			}
		}
	}

  trn.time.start<-proc.time()
	if(cv>1){
		# parallel training to get cv accuracy, # is.ts=TRUE for ts data
		cv.list<-my_cv_partition(trn,is.ts=is.ts,cv=cv)
		
		measure.grid<-array(NA,dim=c(score_grid_length,cv))
		for(i in 1:score_grid_length){
			c=param_c[i]
			d=param_d[i]
			e=param_e[i]
			g=param_g[i]
			for (j in (1:cv.list$k)){ # loop for k fold
				trn.data<-cv.list$trn.data[[j]]
				trn.labels<-cv.list$trn.labels[[j]]
				val.data<-cv.list$val.data[[j]]
				val.labels<-cv.list$val.labels[[j]]
				if(!is.vector(trn.labels)&&!is.factor(trn.labels)){
					trn.labels=trn.labels[,1]
					val.labels=val.labels[,1]
				}
				m<-my_svm(trn.data, trn.labels, x.test=val.data, type=type, kernel=kernel, cost=c, degree=d, epsilon=e, gamma=g, ...)
				pred<-m$predict
				if(type!='C-classification'){
					
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
		best.param.idx<-which(measure.score==min(measure.score),arr.ind=T)
		best.param<-c(param_c[best.param.idx[1]], param_g[best.param.idx[1]], param_e[best.param.idx[1]],param_d[best.param.idx[1]])
		
	}else{# no need cv
		measure.grid<-best.param.idx<-NULL
		best.param<-c(c_range[1],g_range[1],e_range[1],d_range[1])
	}
	names(best.param)<-c('c','g','e','d')
	
	# retrain with best param
	# train initial models
	if(!is.vector(trn$labels)&&!is.factor(trn$labels)){
		trn.labels=trn$labels[,1]
		horizon<-ncol(trn$labels)
	}else{ # 1 ahead forecasting or classification
		trn.labels=trn$labels
		horizon=1
	}
	best.model<-my_svm(trn$data,trn.labels, type=type, kernel=kernel, 
									cost=best.param["c"], degree=best.param["d"], epsilon=best.param["e"], gamma=best.param["g"],...)$model
	trn.time<-proc.time()-trn.time.start
	
	if(is.null(tst)){
		return(list(measure.grid=measure.grid, best.param=best.param, best.param.idx=best.param.idx, model=best.model, predict=NULL,trn.time=trn.time,tst.time=NULL))
	}
  tst.time.start<-proc.time()
	# retrain when new tst data arrived? online training
	tst.data.length<-nrow(tst$data)
	trn.data.length<-nrow(trn$data)
	data.width<-ncol(trn$data)
	
	tst.pred<-array(NA,dim=c(tst.data.length,horizon))
	re_model<-list()
	length(re_model)<-horizon
	counter=1
	
	if(!is.online){# offline training
		if(is.roll){# roll forecast
			tmp<-my.svm.roll.forecast(best.model,trn,tst,is.reest,type,...) # call roll forecast
		}else{#direct forecast
			tmp<-my.svm.direct.forecast(best.model,trn,tst,type,...) # call direct forecast
		}
		re_model<-tmp$model
		tst.pred<-tmp$predict
	}else{# online training
		for (i in 1:tst.data.length){
			# reconstruction tst.data
			trn.data<-rbind(trn$data[i:trn.data.length,],tst$data[(1:i)-1,])
			trn.labels<-rbind(trn$labels[i:trn.data.length,],tst$labels[(1:i)-1,])
			tst.data<-t(as.matrix(tst$data[i,]))
			tst.labels<-t(as.matrix(tst$labels[i,]))
			online.trn<-list(data=trn.data,labels=trn.labels)
			online.tst<-list(data=tst.data,labels=tst.labels)
			
			# check if block update or stepwise update 
			if(is.numeric(is.online)){#block update
				if(counter%%round(n/is.online)==0){ # time to retrain
					if(is.roll){# roll forecast
						re_model<-my.svm.roll.forecast(best.model,online.trn,online.tst,is.reest,type,...)$model
					}else{ # direct forecast
						re_model<-my.svm.direct.forecast(best.model,online.trn,online.tst,type,...)$model
					}
				}else{ # not time to update model
					for (h in 1:horizon){
						re_model[[h]]<-best.model
					}
				}
			}else{# stepwise, always retain
				if(is.roll){# roll forecast
					re_model<-my.svm.roll.forecast(best.model,online.trn,online.tst,is.reest,type,...)$model
				}else{ # direct forecast
					re_model<-my.svm.direct.forecast(best.model,online.trn,online.tst,type,...)$model
				}
			}
			
			for (h in 1:horizon){
				tst.pred[i,h]<-predict(re_model[[h]],online.tst$data)
			}
			counter=counter+1 
		}
	}
	
  tst.time<-proc.time()-tst.time.start
	return(list(measure.grid=measure.grid, best.param=best.param, best.param.idx=best.param.idx, model=re_model, predict=tst.pred, trn.time=trn.time,tst.time=tst.time))
}




# build in functions for my_cv_svm
my.svm.roll.forecast<-function(model,trn,tst,is.reest=FALSE,type,...){
	if(model$kernel==0){
		kernel='linear'
	}else if(model$kernel==1){
		kernel='polynomial'
	}else if(model$kernel==2){
		kernel='radial'
	}else{
		kernel='sigmoid'
	}
	
	# rolling forecast
	horizon<-ncol(trn$labels)
	# initialize
	h.model<-list()
	length(h.model)<-horizon
	h.fit<-array(0,dim=dim(trn$labels))
	h.pred<-array(0,dim=dim(tst$labels))
	
	# loop for each horizion
	for (h in 1:horizon){
		h.trn.labels<-trn$labels[,h]# 1 step ahead first
		h.tst.labels<-tst$labels[,h]# 1 step ahead first
		if(h==1){ # first time
			#data.width<-ncol(trn$data)
			h.trn.data<-trn$data
			h.tst.data<-tst$data
			h.model[[h]]=model
		}else{ # h>1
			# update new dataset for training
			h.trn.data<-cbind(matrix(h.trn.data[,-1],nrow=nrow(trn$data)),h.fit[,h-1])
			h.tst.data<-cbind(matrix(h.tst.data[,-1],nrow=nrow(tst$data)),h.pred[,h-1])
			if(is.reest){#update on each roll
				h.model[[h]]=my_svm(h.trn.data, h.trn.labels, type=type, kernel=kernel, 
												 cost=model$cost, degree=model$degree, epsilon=model$epsilon, gamma=model$gamma,...)$model
			}else # keep using same model
			{
				h.model[[h]]=model
			}
		}
		
		#browser()
		# forecast
		h.pred[,h]<-predict(h.model[[h]],h.tst.data)
		h.fit[,h]<-predict(h.model[[h]],h.trn.data)
	}
	return(list(model=h.model,predict=h.pred,fit=h.fit))
}

my.svm.direct.forecast<-function(model,trn,tst,type,...){
	if(model$kernel==0){
		kernel='linear'
	}else if(model$kernel==1){
		kernel='polynomial'
	}else if(model$kernel==2){
		kernel='radial'
	}else{
		kernel='sigmoid'
	}
	
	# direct forecast
	# initialize
	h.model<-list()
	horizon<-ncol(as.matrix(trn$labels))
	if(horizon==1){#only one horizon
		h.fit<-rep(0,length(trn$labels))
		h.pred<-rep(0,length(tst$labels))
	}else{
		h.fit<-array(0,dim=dim(trn$labels))
		h.pred<-array(0,dim=dim(tst$labels))
	}
	
	# loop for each horizion
	for (h in 1:horizon){
		if(horizon==1)
			h.trn.labels<-trn$labels
		else
			h.trn.labels<-trn$labels[,h]# 1 step ahead first
		
		h.model[[h]]=my_svm(trn$data, h.trn.labels, type=type, kernel=kernel, 
										 cost=model$cost, degree=model$degree, epsilon=model$epsilon, gamma=model$gamma,...)$model
		
		
		# forecast
		tmp.pred<-predict(h.model[[h]],tst$data)
		tmp.fit<-predict(h.model[[h]],trn$data)
		if(horizon==1){ # no need roll, only one data
			h.pred<-tmp.pred
			h.fit<-tmp.fit
		}else{ # multiple horizon
			h.pred[,h]<-tmp.pred
			h.fit[,h]<-tmp.fit
		}
	}
	return(list(model=h.model,predict=h.pred,fit=h.fit))
}















# bagging svr
my_bagging_svr<-function(trn,tst,LAG=2,B=100,SD=0.2,
												 kernel='radial',
												 c_range=2^(-8:8),
												 d_range=2:5,
												 g_range=(1/LAG)*(1:(LAG-1)),
												 e_range=2^(-8:0),k=5,is.ts=FALSE) {
	
	#   require(e1071)
	#   require(bootstrap)
	#   require(foreach)
	#   require(doParallel)
	
	#for ts data, need to do cv_ts partition first
	if(is.ts){
		ts_trn<-my_ts_cv_partition<-function(trn,cv=5)
			
			# coarse grid search
			best.param.list<-my_svr_cv.parallel(ts_trn,c_range=c_range,g_range=g_range,e_range=e_range,k=k)
		best.param<-best.param.list$best.param
		m.coarse<-svm(trn$data,trn$labels,cost=best.param[1],gamma=best.param[2],epsilon=best.param[3])
		# refine SVR parameter selection
		c_range_fine<-2^(log2(best.param[1])+(-2:2))
		g_range_fine<-best.param[2]*seq(0.2:1, by=0.2)
		e_range_fine<-2^((log2(best.param[3]))+(-4:0))
		d_range_fine<-best.param[4]
		# fine grid search
		best.param.fine.list<-my_svr_cv.parallel(ts_trn,kernel='radial',c_range_fine,g_range_fine,e_range_fine,d_range_fine, k)
		best.param.fine<-best.param.fine.list$best.param
		m.fine<-svm(trn$data, trn$labels,cost=best.param.fine[1],gamma=best.param.fine[2],epsilon=best.param.fine[3])
		
	}else
	{
		
		# coarse grid search
		best.param.list<-my_svr_cv.parallel(x=trn,c_range=c_range,g_range=g_range,e_range=e_range,k=k)
		best.param<-best.param.list$best.param
		m.coarse<-svm(trn$data,trn$labels,cost=best.param[1],gamma=best.param[2],epsilon=best.param[3])
		# refine SVR parameter selection
		c_range_fine<-2^(log2(best.param[1])+(-2:2))
		g_range_fine<-best.param[2]*seq(0.2:1, by=0.2)
		e_range_fine<-2^((log2(best.param[3]))+(-4:0))
		d_range_fine<-best.param[4]
		# fine grid search
		best.param.fine.list<-my_svr_cv.parallel(trn,kernel='radial',c_range_fine,g_range_fine,e_range_fine,d_range_fine, k)
		best.param.fine<-best.param.fine.list$best.param
		m.fine<-svm(trn$data, trn$labels,cost=best.param.fine[1],gamma=best.param.fine[2],epsilon=best.param.fine[3])
	}
	
	# bootstrap
	idx.length<-nrow(trn$data)
	boot<-bootstrap(1:idx.length, B, as.numeric)
	idx.boot<-boot$thetastar # oob.idx=-idx.boot
	
	# create randomness to the best parameters
	rand.param<-rnorm(4*B, mean=rep(best.param.fine, each=B), sd=SD)
	rand.param<-array(rand.param, dim=c(B,4))
	rand.param<-abs(rand.param) # delete negative 
	rand.param[which(rand.param[,3]>1),3]<-1 # epsilon less than 1
	colnames(rand.param)<-c('c','g','e','d')
	
	model.list<-list()
	length(model.list)<-B
	single.pred<-array(0, dim=c(nrow(tst$data), B))
	single.error<-array(0, dim=c(B,3))
	oob.error<-array(0, dim=c(B,3))
	# for each bag
	for (j in (1:B)) # change to parallel
	{
		trn.bag<-trn$data[idx.boot[,j],] # training data
		trn.oob<-trn$data[-idx.boot[,j],] # oob data
		trn.bag.labels<-trn$labels[idx.boot[,j]] # trining labels
		trn.oob.labels<-trn$labels[-idx.boot[,j]] # oob labels
		m<-svm(trn.bag, trn.bag.labels, kernel='radial', type='eps-regression', cost=rand.param[j,1], epsilon=rand.param[j,3])
		model.list[[j]]<-m
		oob.pred<-predict(m, trn.oob)  
		
		oob<-my_error_measure(oob.pred, trn.oob.labels, trn.oob[,ncol(trn.oob)])
		oob.error[j,]<-c(oob$RMSE, oob$MAPE, oob$MASE)
		single.pred[,j]<-predict(m, tst$data)
		single<-my_error_measure(single.pred[,j], tst$labels,tst$data[,ncol(tst$data)])
		single.error[j,]<-c(single$RMSE, single$MAPE, single$MASE)
	}
	colnames(oob.error)<-c('RMSE','MAPE','MASE')
	colnames(single.error)<-c('RMSE','MAPE','MASE')
	
	# compare bagging and single
	bagging.pred<-rowMeans(single.pred)
	coarse.pred<-predict(m.coarse, tst$data)
	fine.pred<-predict(m.fine, tst$data)
	
	bagging.error<-my_error_measure(bagging.pred, tst$labels, err.benchmk=tst$data[,ncol(tst$data)])
	coarse.error<-my_error_measure(coarse.pred, tst$labels, err.benchmk=tst$data[,ncol(tst$data)])
	fine.error<-my_error_measure(fine.pred, tst$labels, err.benchmk=tst$data[,ncol(tst$data)])
	
	return(list(bagging.model=model.list, coarse.model=m.coarse, fine.model=m.fine,
							best.param.fine.list=best.param.fine.list,
							best.param.coarse.list=best.param.list,
							idx.boot=idx.boot, single.pred=single.pred, bagging.pred=bagging.pred,
							coarse.pred=coarse.pred, fine.pred=fine.pred, oob.error=oob.error,
							single.error=single.error, bagging.error=bagging.error, coarse.error=coarse.error,
							fine.error=fine.error))
}
AdaBoost.M1<-function( x, y, x.test=NULL, method='mlp', maxIter=100, ... ){
	# Experiments with a New Boosting Algorithm
	# Yoav Freund Robert E. Schapire (1996), ICML
	
	# example
	# 	data("iris")
	# 	iris <- iris[sample(1:nrow(iris), length(1:nrow(iris))), 1:ncol(iris)]
	# 	irisValues <- iris[, 1:4]
	# 	irisTargets <- iris[, 5]
	# 	iris <- splitForTrainingAndTest(irisValues, irisTargets, ratio = 0.15)
	# 	iris <- normTrainingAndTestSet(iris)
	# 	AdaBoost.M1(iris$inputsTrain, iris$targetsTrain, method='mlp', maxIter=100, x.test = iris$inputsTest,  size = 9,	learnFuncParams = 0.1, maxit = 60 )
	
	
	y<-decodeClassLabels(y) # convert to nominal type
	class.labels<-colnames(y) # categories of y
	k<-ncol(y) # number of classes
	m<-nrow(x) # parameter of x
	n<-ncol(x)
	#initialize weight distr and error
	D_new=rep(1/m, m)
	
	T=maxIter # max iter
	eth=0.5 # error threshold
	
	fit<-D<-array(NA,dim=c(T,m))
	e=rep(0,T)
	beta=rep(1,T)
	model=list()
	# iterate
	t=1
	th_e=0
	runtime.counter=1
	while (t<=T && runtime.counter<100){
		D[t,]=D_new
		# bootstrap with dist D
		idx<-sample(m,replace=TRUE,prob=D_new)
		tth_x=x[idx,]
		tth_y=y[idx,]
		if(ncol(unique(tth_y))==1){ # single class, redo
			runtime.counter=runtime.counter+1
			next
		}
		# call weak learner, CART
		model[[t]]=eval(parse(text=method))(tth_x,tth_y, ...)
		# getback to hypothesis
		fit[t,]<-y_fit<-encodeClassLabels(predict(model[[t]],x))
		# calculate error rate
		th_e<-e[t]<-sum(D_new%*%as.numeric(y_fit!=encodeClassLabels(y)))
		if(e[t]==0) # perfect, bypass this one, hardly happen
			next
		else if(e[t]>eth) # more than threshold 0.5 by default
			break
		
		# question!!!
		# if all correct (perfect predict), then e=0 and beta=0
		
		# weight updating param
		beta[t]=th_e/(1-th_e)
		# update distr
		correct=beta[t]*as.numeric(y_fit==encodeClassLabels(y))
		incorrect=as.numeric(y_fit!=encodeClassLabels(y))
		D_new=D_new*(correct+incorrect)
		# normalize
		D_new=D_new/sum(D_new)
		t=t+1
	}
	max_t=t-1
	
	# final model: cart, beta(t)
	ww=log(1/beta[1:max_t])
	
	# final fit value
	aggregated.fit<-my_majority_vote(fit[1:max_t,],weight=ww,levels=1:k,labels=class.labels,byrow=FALSE)
	
	
	# predict
	if(!is.null(x.test)){ # have test dataset
		m.test<-nrow(x.test)
		predict<-array(NA,dim=c(max_t,m.test))
		for (t in 1:max_t){
			predict[t,]<-encodeClassLabels(predict(model[[t]],x.test))
		}
		aggregated.predict<-my_majority_vote(predict,weight=ww,levels=1:k,labels=class.labels,byrow=FALSE)
	}else{
		aggregated.predict=NULL
	}
	
	return(list(fit=aggregated.fit$value,predict=aggregated.predict$value,beta=beta[1:max_t],weight=ww,max.t=max_t))
}
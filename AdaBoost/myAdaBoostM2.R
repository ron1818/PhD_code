# incomplete


AdaBoost.M2<-function( x, y, x.test=NULL, method='mlp', maxIter=100,... ){
	# Analysis of the Performance of AdaBoost.M2 for
	# the Simulated Digit-Recognition-Example
	# G¨unther Eibl and Karl Peter Pfeiffer
	
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
	G<-ncol(y) # number of classes
	N<-nrow(x) # parameter of x
	n<-ncol(x)
	# initialize y_fit, all incorrect
	y_fit<-matrix(as.numeric(!y),ncol=G)
	# initialize weight distr and error
	D_new=rep(1/N, N)
	w_new=D_new/(G-1) # weight vector
	T=maxIter # max iter
	
	W<-w<-fit<-D<-array(NA,dim=c(T+1,N))
	e=rep(0,T)
	alpha=rep(1,T)
	model=list()
	# iterate
	t=1
	th_e=0
	runtime.counter=1
	q<-rep(1,N)
	while (t<=T){
		w[t,]=w_new
		incorrect.sample=encodeClassLabels(y_fit)!=encodeClassLabels(y)
		W[t,]=as.numeric(incorrect.sample)*w_new
		q[incorrect.sample]=w[t,incorrect.sample]/W[t,incorrect.sample]
		D[t,]=W[t,]/sum(W[t,])
		
		# bootstrap with dist D
		idx<-sample(m,replace=TRUE,prob=D[t,])
		tth_x=x[idx,]
		tth_y=y[idx,]
		if(ncol(unique(tth_y))==1){ # single class, redo
			runtime.counter=runtime.counter+1
			next
		}
		# call weak learner, CART
		model[[t]]=eval(parse(text=method))(tth_x,tth_y)
		# getback to hypothesis, probability
		y_fit<-predict(model[[t]],x)
		incorrect.sample=encodeClassLabels(y_fit)!=encodeClassLabels(y)
		# calculate pseudoloss
		correct.prob<-apply(y_fit*y,1,sum) # ht(xi,yi)
		incorrect.prob<-apply(y_fit*!y,1,sum) # ht(xi,y)
		th_e<-e[t]<-0.5*sum(D[t,]*(1-correct.prob+1/(G-1)*as.numeric(incorrect.sample)*q*(incorrect.prob)))
		alpha[t]=0.5*log((1-th_e)/th_e)
		# update distr
		w[t+1,incorrect.sample]<-w[t,incorrect.sample]*exp(-alpha[t]*(1-correct.prob[incorrect.sample]+incorrect.prob[incorrect.sample]))
		t=t+1
	}

}
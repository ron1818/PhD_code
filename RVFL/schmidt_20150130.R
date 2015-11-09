

# define function sigmoid
sigmoid<-function(x){
	x<-as.matrix(x)
	# output is from 0 to 1
	return(1/(1+exp(-x)))
}

#Schmidt1992
# random weights on hidden layer M, one bias node
# least square on output layer weights W, one bias node
# no direct link from input to output
# input x, output y
# input to hidden weights: M
# use runif -1 1 for M
# pseudo inverse for W
my.schmidt.rvfl<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL){
  # dim information
	N<-nrow(x)
	n_i<-ncol(x)
  # initialize test parameter 
  mse.fit<-rep(Inf, length(n_h))
	if(!is.null(x.test)){
    mse.pred<-rep(Inf, length(n_h))
  	y.pred<-array(NA,dim=c(nrow(x.test),length(n_h)))
	}else{
    mse.pred<-NULL
    y.pred<-NULL
	}
  # initialize train parameter  
	W<-list()
	length(W)<-length(n_h)
  # create randome weights
	tmp<-runif(max(n_h)*(n_i+1), -1, 1)
	M<-matrix(tmp,nrow=max(n_h),ncol=n_i+1)
	c=1
	for (i in n_h){
# 		# M, random generated
# 		tmp<-runif(i*(n_i+1), -1, 1)
# 		M[[c]]<-matrix(tmp,nrow=i,ncol=n_i+1)
		
		# activation function
		x.with.bias<-cbind(x,rep(1,N))
		o.hidden<-sigmoid(M[1:i,]%*%t(x.with.bias))
		o.hidden.with.bias<-t(rbind(o.hidden,rep(1,N)))
		W[[c]]<-ginv(o.hidden.with.bias)%*%y
    mse.fit[c]<-mean((y-o.hidden.with.bias%*%W[[c]])^2)
    
		### testing
		if(!is.null(x.test)){
      N.t<-nrow(x.test)
		  x.test.with.bias<-cbind(x.test,rep(1,N.t))
      
		  o.hidden.t<-sigmoid(M[1:i,]%*%t(x.test.with.bias))
		  o.hidden.t.with.bias<-t(rbind(o.hidden.t,rep(1,N.t)))
      y.pred[,c]<-o.hidden.t.with.bias%*%W[[c]]
		  mse.pred[c]<-mean((y.test-y.pred[,c])^2)
		}
    
		c=c+1
	}
	return(list(W=W,M=M,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred))
}

#Pao 1994
# bias on each enhancement node
# direct link from input to output
my.pao.rvfl<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL){
  # dim information	
	N<-nrow(x)
	n_i<-ncol(x)
	# initialize test parameter  
	mse.fit<-rep(Inf, length(n_h))
	if(!is.null(x.test)){
	  mse.pred<-rep(Inf, length(n_h))
	  y.pred<-array(NA,dim=c(nrow(x.test),length(n_h)))
	}else{
	  mse.pred<-NULL
	  y.pred<-NULL
	}
	# initialize train parameter    
	# create randome weights
	tmp<-runif(max(n_h)*(n_i), -1, 1)
	A<-matrix(tmp,nrow=max(n_h),ncol=n_i)
	W<-list()
	length(W)<-length(n_h)
	c=1
	for (i in n_h){
# 		# A, random generated
# 		tmp<-runif(i*n_i, -1, 1)
# 		A[[c]]<-matrix(tmp,nrow=i,ncol=n_i)
# 		# b, threshold parameter
    a.times.x<-A[1:i,]%*%t(x)
		b<--1*apply(a.times.x,1,mean)
		# enhancement node
		enhancement.node<-t(sigmoid(a.times.x+b))
		# all nodes
		all.node<-cbind(x,enhancement.node)
		W[[c]]<-ginv(all.node)%*%y
    #mse.fit[c]<-mean((y-all.node%*%W[[c]])^2)
    
		### testing
		if(!is.null(x.test)){
		  a.times.x.t<-A[1:i,]%*%t(x.test)
		  # enhancement node
		  en.t<-t(sigmoid(a.times.x.t+b))
		  # all nodes
		  all.node.t<-cbind(x.test,en.t)
		  y.pred[,c]<-all.node.t%*%W[[c]]
		  mse.pred[c]<-mean((y.test-y.pred[,c])^2)
		}
		c=c+1
	}
	return(list(W=W,A=A,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred))
}

#Chen 1999
# bias on each enhancement node
# direct link from input to output
my.chen.rvfl<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL){

  # dim information  
  N<-nrow(x)
  n_i<-ncol(x)
  # initialize test parameter  
  mse.fit<-rep(Inf, length(n_h))
  if(!is.null(x.test)){
    mse.pred<-rep(Inf, length(n_h))
    y.pred<-array(NA,dim=c(nrow(x.test),length(n_h)))
  }else{
    mse.pred<-NULL
    y.pred<-NULL
  }
  # initialize train parameter    
  # create randome weights
  tmp<-runif(max(n_h)*(n_i), -1, 1)
  A<-matrix(tmp,nrow=max(n_h),ncol=n_i)
  W<-list()
  length(W)<-length(n_h)
  c=1
  for (i in n_h){
    ## A, random generated
    # tmp<-runif(i*n_i, -1, 1)
    # [[c]]<-matrix(tmp,nrow=i,ncol=n_i)
    # b, threshold parameter
    a.times.x<-A[1:i,]%*%t(x)
    b<-runif(i,-1,1)
    # enhancement node
    enhancement.node<-t(tanh(a.times.x+b))
    # all nodes
    all.node<-cbind(x,enhancement.node)
    W[[c]]<-ginv(all.node)%*%y
    #mse.fit[c]<-mean((y-all.node%*%W[[c]])^2)
    
    ### testing
    if(!is.null(x.test)){
      a.times.x.t<-A[1:i,]%*%t(x.test)
      # enhancement node
      en.t<-t(tanh(a.times.x.t+b))
      # all nodes
      all.node.t<-cbind(x.test,en.t)
      y.pred[,c]<-all.node.t%*%W[[c]]
      mse.pred[c]<-mean((y.test-y.pred[,c])^2)
    }
    c=c+1
  }
  return(list(W=W,A=A,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred))
}


# k-fold CV version of the RVFLs
my.cv.pao.rvfl<-function(x,y,n_h_range,cv=5,x.test=NULL,y.test=NULL){
  # cross validation
  # cv to find optimal # of enhancement node
  trn<-list(data=x,labels=y) # put into list
  cv.list<-my_ts_cv_partition(trn,cv=cv)
  cv.mse.grid<-foreach(c=1:cv, .combine='rbind',.export=c('sigmoid','my.pao.rvfl'),.packages='MASS')%dopar%{
    x<-cv.list$trn.data[[c]]
    y<-cv.list$trn.labels[[c]]
    x.test<-cv.list$val.data[[c]]
    y.test<-cv.list$val.labels[[c]]
    
    cv.model<-my.pao.rvfl(x,y,n_h_range,x.test,y.test)
    return(cv.model$mse.pred)
  }
  cv.mse.grid<-apply(cv.mse.grid,2,mean)
  best.n_h<-n_h_range[which(cv.mse.grid==min(cv.mse.grid,na.rm=T))]
  
  best.model<- my.pao.rvfl(x,y,best.n_h,x.test,y.test)
  predict=best.model$predict[[1]]
  return(list(model=best.model,predict=predict,mse.grid=cv.mse.grid,cv=cv))
}

my.cv.chen.rvfl<-function(x,y,n_h_range,cv=5,x.test=NULL,y.test=NULL){
    # change x from 0 1 to -1 1
  x<-2*x-1
  y<-2*y-1
  if(!is.null(x.test)){
    x.test<-2*x.test-1
    y.test<-2*y.test-1
  }
  
  # cross validation
  # cv to find optimal # of enhancement node
  trn<-list(data=x,labels=y) # put into list
  cv.list<-my_ts_cv_partition(trn,cv=cv)
  cv.mse.grid<-foreach(c=1:cv, .combine='rbind',.export=c('sigmoid','my.chen.rvfl'),.packages='MASS')%dopar%{
    x<-cv.list$trn.data[[c]]
    y<-cv.list$trn.labels[[c]]
    x.test<-cv.list$val.data[[c]]
    y.test<-cv.list$val.labels[[c]]
    
    cv.model<-my.chen.rvfl(x,y,n_h_range,x.test,y.test)
    return(cv.model$mse.pred)
  }
  cv.mse.grid<-apply(cv.mse.grid,2,mean)
  best.n_h<-n_h_range[which(cv.mse.grid==min(cv.mse.grid,na.rm=T))]
  
  best.model<- my.chen.rvfl(x,y,best.n_h,x.test,y.test)
  predict=best.model$predict[[1]]
  return(list(model=best.model,predict=predict,mse.grid=cv.mse.grid,cv=cv))
}

my.cv.schmidt.rvfl<-function(x,y,n_h_range,cv=5,x.test=NULL,y.test=NULL){
  # cross validation
  # cv to find optimal # of enhancement node
  trn<-list(data=x,labels=y) # put into list
  cv.list<-my_ts_cv_partition(trn,cv=cv)
  cv.mse.grid<-foreach(c=1:cv, .combine='rbind',.export=c('sigmoid','my.schmidt.rvfl'),.packages='MASS')%dopar%{
    x<-cv.list$trn.data[[c]]
    y<-cv.list$trn.labels[[c]]
    x.test<-cv.list$val.data[[c]]
    y.test<-cv.list$val.labels[[c]]
    
    cv.model<-my.schmidt.rvfl(x,y,n_h_range,x.test,y.test)
    return(cv.model$mse.pred)
  }
  cv.mse.grid<-apply(cv.mse.grid,2,mean)
  best.n_h<-n_h_range[which(cv.mse.grid==min(cv.mse.grid,na.rm=T))]
  
  best.model<- my.schmidt.rvfl(x,y,best.n_h,x.test,y.test)
  predict=best.model$predict[[1]]
  return(list(model=best.model,predict=predict,mse.grid=cv.mse.grid,cv=cv))
}
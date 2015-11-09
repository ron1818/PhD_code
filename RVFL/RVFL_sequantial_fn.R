my.RVFL.ts.addnode<-function(A,a,Y,G=NULL){
  if(!is.null(G)){#weighted LSE
    #    browser()
    A<-sqrt(G)%*%A
    Y<-sqrt(G)%*%Y
  }
  
  # pseudoinverse of A_n
  pseudoinv.A<-ginv(A)
  # W<-qr.solve(A,Y)
  W<-pseudoinv.A%*%Y
  # calculate d
  d<-pseudoinv.A%*%a
  # calculate c
  c=a-A%*%d
  # calculate b'
  if(all(c<=1e-8)){
    b=solve(1+t(d)%*%d)%*%t(d)%*%pseudoinv.A
  }else{
    b=ginv(c)
  }
  W_new<-rbind(W-d%*%b%*%Y,b%*%Y)
  pseudoinv.A_new<-rbind(pseudoinv.A-d%*%b,b)
  A_new=ginv(pseudoinv.A_new)
  return(list(A=A_new,W=W_new))
}

my.RVFL.ts.addsample<-function(A,a,Y,y, G=NULL){
  if(!is.null(G)){#weighted LSE
    #    browser()
    A<-sqrt(G)%*%A
    Y<-sqrt(G)%*%Y
  }
  # pseudoinverse of A_n
  pseudoinv.A<-ginv(A)
  # W<-qr.solve(A,Y)
  W<-pseudoinv.A%*%as.matrix(Y)
  A_new<-rbind(A,a)
  # calculate d'
  # browser()
  d.prime<-a%*%pseudoinv.A
  # calculate c
  c=a-d.prime%*%A
  # calculate b'
  if(all(c<1e-8)){
    b=(as.numeric(1+d.prime%*%t(d.prime)))^(-1)*pseudoinv.A%*%t(d.prime)
  }else{
    b=ginv(t(c))
  }
  pseudoinv.A_new<-cbind(pseudoinv.A-b%*%d.prime,b)
  W_new<-W+b%*%(t(y)-t(a)%*%W)
  return(list(A=A_new,W=W_new))
}

my.RVFL.ts.weights<-function(Y,Y_pred,rho=0.5){
  # weighted version
  R<-Y-Y_pred # residual
  M<-median(abs(R))
  rho=0.5 # control shape of RBF
  G<-exp(-(R/(rho*M))^2) # RBF function
  G<-diag(as.numeric(G)) # diagnoal matrix
  return(G)
}

my.RVFL.stepwise<-function(x,y,enhancement.node=nrow(x)-ncol(x),mse.stop=0.01,stepwise=TRUE,x.test=NULL,y.test=NULL)
{
  # train RVFL
  # for stepwise add enhancement node
  max.iter=enhancement.node
  # some parameters
  N<-nrow(x)
  # add bias to training data
  x<-cbind(x,rep(1,N))
  n<-ncol(x)
  # enhancement node's random weights
  w_h<-runif(n*(max.iter),min(x)/n,max(x)/n) # batch create w_h, avoid saturation
  w_h<-matrix(w_h,nrow=n,ncol=max.iter)
  # initialize
  mse.fit<-rep(Inf,max.iter+1)
  A<-W<-list()
  length(A)<-length(W)<-max.iter+1
  A[[1]]<-x # first is the one w/o enhancement
  W[[1]]<-ginv(A[[1]])%*%y # first is the one w/o enhancement
  mse.fit[1]<-mean((y-A[[1]]%*%W[[1]])^2)
  
  if(!is.null(x.test)){ # if have test data
    N.test<-nrow(x.test)
    x.test<-cbind(x.test,rep(1,N.test)) # add bias to test data
    mse.val<-rep(Inf,max.iter+1)
    mse.val[1]<-mean((y.test-x.test%*%W[[1]])^2)
  }else{
    mse.val=NULL
  }
  if (!stepwise){# non stepwise
    a<-tanh(A[[1]]%*%w_h)
    A[[2]]<-cbind(A[[1]],a)
    W[[2]]<-ginv(A[[2]])%*%y
    mse.fit[2]<-mean((y-A[[2]]%*%W[[2]])^2) # calculate mse
    if (!is.null(x.test)){ # if have test data
      test.new<-cbind(x.test,tanh(x.test%*%w_h))
      mse.val[2]<-mean((y.test-test.new%*%W[[2]])^2)
    }
  }else{ # stepwise
    # initialize
    iter=1
    terminate.counter=1
    TERMINATE=5
    current.A<-x
    current.mse=Inf
    A_added<-NULL
    A_v_added<-NULL
    while(current.mse>mse.stop & iter<=max.iter){ # loop until stop criteria
      a<-tanh(A[[1]]%*%w_h[,iter]) # add one hidden node
      
      if(qr(current.A)$rank==qr(cbind(current.A,a))$rank){# rank not increased
        if(terminate.counter>=TERMINATE){
          break
        }else{
          terminate.counter=terminate.counter+1
          w_h[,iter]<-runif(n,min(x)/n,max(x)/n)
          w_h[,iter]<-matrix(w_h[,iter],nrow=n,ncol=1)
          next # update a new w_h and do it again 
        }
      }else{
        terminate.counter=1 # reset terminate.counter
        RVFL.m<-my.RVFL.ts.addnode(current.A,a,y,G=NULL)
        A[[iter+1]]<-current.A<-RVFL.m$A # update A
        W[[iter+1]]<-current.W<-RVFL.m$W # update W
        current.w_h<-w_h[,iter]
      }
      
      # fit prediction
      fit<-current.A%*%current.W
      mse.fit[iter+1]<-current.mse<-mean((y-fit)^2)
      
      if (!is.null(x.test)){# with test data
        # update val data with enhancement node
        av<-tanh(x.test%*%current.w_h) # add one hidden nodes
        A_v_added<-cbind(A_v_added,av)
        A_v_new<-cbind(x.test,A_v_added)
        # prediction
        #browser()
        pred<-A_v_new%*%current.W
        mse.val[iter+1]<-mean((y.test-pred)^2)
      }
      iter=iter+1
#  browser()
    }
  }  
  return(list(A=A,W=W,w_h=w_h,mse.val=mse.val,mse.fit=mse.fit))
}

my.RVFL.cv<-function(x,y,enhancement.node=nrow(x)-ncol(x),mse.stop=0.01,stepwise=TRUE,cv=5,x.test=NULL,y.test=NULL){
  # cross validation
  # cv to find optimal # of enhancement node
  trn<-list(data=x,labels=y) # put into list
  cv.list<-my_ts_cv_partition(trn,cv=cv)
  cv.mse.grid<-foreach(c=1:cv, .combine='rbind',.export=c('my.RVFL','my.RVFL.ts.addnode','my.RVFL.ts.addsample'),.packages='MASS')%dopar%{
    x<-cv.list$trn.data[[c]]
    y<-cv.list$trn.labels[[c]]
    x.test<-cv.list$val.data[[c]]
    y.test<-cv.list$val.labels[[c]]
    
    cv.model<-my.RVFL(x,y,x.test=x.test,y.test=y.test)
    mse<-cv.model$mse.val
    return(mse)
  }
  cv.mse.grid<-apply(cv.mse.grid,2,mean)
  best.enhancement.node<-which(cv.mse.grid==min(cv.mse.grid,na.rm=T))-1
  
  if(!is.null(x.test)){
    # use best.idx to train on whole trn dataset
    best.model<-my.RVFL(x,y,enhancement.node=best.enhancement.node,stepwise=F)
    x.test<-cbind(x.test,rep(1,nrow(x.test)))
    add.test<-tanh(x.test%*%best.model$w_h)
    test.pred<-cbind(x.test,add.test)%*%best.model$W[[2]]
  }else{
    test.pred=NULL
  }
  return(list(model=best.model,predict=test.pred,best.enhancement.node=best.enhancement.node))
}
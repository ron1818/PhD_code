# Chen, RVFL, IEEE TNN, 1996
library(MASS)
library(forecast)

source('pre_processing_fn.R')
source('post_processing_fn.R')
source('RVFL//RVFL_fn.R')
source('cross_validation_fn.R')

org.ts<-arima.sim(n=500,list(ar=c(0.88,-0.48)),sd=sqrt(0.18))
ts<-my_scale_ts(org.ts,c(0,1))
tmp<-my_vector_split(ts$scaled.ts,HORIZON=1,LAG=2,RATIO=0.5,point=F)
train<-tmp$trn
test<-tmp$tst

x<-train$data
x.test<-test$data
y<-train$labels
y.test<-test$labels


A_train<-cbind(train$data,rep(1,nrow(train$data))) # add bias
A_test<-cbind(test$data,rep(1,nrow(test$data))) # add bias
y_train<-train$labels
y_test<-test$labels

# some parameters
N_n<-nrow(A_train)
n_n<-ncol(A_train) # with bias
r_n<-qr(A_train)$rank
m_n<-ncol(y_train)


# # Initialize RVFL without enhancement node
# W_init<-qr.solve(A_train,y_train)

# CV Train recursively with adding enhancement node
cv=5
cv.data<-my_ts_cv_partition(list(data=A_train,labels=y_train),cv=cv)
mse.stop=0.01
max.iter=max(100,nrow(cv.data$trn.data))

# initialize

A<-W<-list()
length(A)<-length(W)<-cv
mse.fit<-array(Inf,dim=c(cv,max.iter+1))
mse.val<-array(Inf,dim=c(cv,max.iter+1))
w_h<-runif(n_n*(max.iter+1),-1,1) # batch create w_h
w_h<-matrix(w_h,nrow=n_n,ncol=max.iter+1)

# cv to find best # enhancement node
for( i in 1:cv ){ # for each fold
  # initialize
  iter=2
  current.mse=Inf
  A_added<-NULL
  A_v_added<-NULL
  
  A[[i]]<-W[[i]]<-list()
  length(A[[i]])<-length(W[[i]])<-max.iter+1
  
  Y<-as.matrix(cv.data$trn.labels[[i]])
  Y.val<-as.matrix(cv.data$val.labels[[i]])
  A[[i]][[1]]<-current.A<-cv.data$trn.data[[i]] # first is the one w/o enhancement
  W[[i]][[1]]<-ginv(A[[i]][[1]])%*%Y # first is the one w/o enhancement
  mse.fit[i,1]<-mean((Y-A[[i]][[1]]%*%W[[i]][[1]])^2)
  mse.val[i,1]<-mean((Y.val-cv.data$val.data[[i]]%*%W[[i]][[1]])^2)

  while(current.mse>mse.stop & iter<=max.iter){ # loop until stop criteria
    a<-tanh(A[[i]][[1]]%*%w_h[,iter]) # add one hidden node
    
    if(qr(current.A)$rank==qr(cbind(current.A,a))$rank){# rank not increased
      w_h[,iter]<-runif(n_n,-1,1)
      w_h[,iter]<-matrix(w_h[,iter],nrow=n_n,ncol=1)
      next # update a new w_h and do it again 
    }else{
      RVFL.m<-my.RVFL.ts.addnode(current.A,a,Y,G=NULL)
      A[[i]][[iter]]<-current.A<-RVFL.m$A # update A
      W[[i]][[iter]]<-current.W<-RVFL.m$W # update W
      current.w_h<-w_h[,iter]
    }
    
    # prediction
    fit<-current.A%*%current.W
    mse.fit[i,iter]<-current.mse<-mean((cv.data$trn.labels[[i]]-fit)^2)
    
    # update val data with enhancement node
    av<-tanh(cv.data$val.data[[i]]%*%current.w_h) # add one hidden nodes
    A_v_added<-cbind(A_v_added,av)
    A_v_new<-cbind(cv.data$val.data[[i]],A_v_added)
    # prediction
    #browser()
    pred<-A_v_new%*%current.W
    mse.val[i,iter]<-mean((cv.data$val.labels[[i]]-pred)^2)
    iter=iter+1
  }
  
}

cv.mse.val<-apply(mse.val,2,mean)
best.idx<-which(cv.mse.val==min(cv.mse.val))

# use best.idx to train on whole trn dataset
best.w_h<-w_h[,1:best.idx]
add.train<-tanh(A_train%*%best.w_h)
best.A<-cbind(A_train,add.train)
best.W<-ginv(best.A)%*%y_train

# predict on test dataset, offline testing
add.test<-tanh(A_test%*%best.w_h)
pred.offline<-cbind(A_test,add.test)%*%best.W


# online testing
pred.online<-rep(NA,length(y_test))
# first sample
a<-A_test[1,]
y<-y_test[1]
pred.online[1]<-c(a,tanh(a%*%best.w_h))%*%best.W

A=best.A
Y=y_train
for (i in 2:length(y_test)){
  a<-A_test[i,]
  a<-c(a,tanh(a%*%best.w_h)) # append with enhancement node
  y<-y_test[i]
  tmp<-my.RVFL.ts.addsample(A,a,Y,y)
  #A<-tmp$A
  W<-tmp$W
  pred.online[i]<-a%*%W
}

offline.error<-my_error_measure(pred.offline,y_test,A_test[,ncol(A_test)-1])
online.error<-my_error_measure(as.matrix(pred.online),y_test,A_test[,ncol(A_test)-1])

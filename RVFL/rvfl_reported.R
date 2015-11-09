
#Schmidt1992
# random weights on hidden layer M, one bias node
# least square on output layer weights W, one bias node
# no direct link from input to output
# input x, output y
# input to hidden weights: M
# use runif -1 1 for M
# pseudo inverse for W
my.schmidt.rvfl<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL,seed=NULL){
    y<-as.matrix(y)
    # dim information
    N<-nrow(x)
    n_i<-ncol(x)
    n_o<-ncol(y)
    mse.fit<-array(Inf, dim=c(length(n_h),n_o))
    # initialize test parameter
    if(!is.null(x.test)){
        y.test<-as.matrix(y.test)
        mse.pred<-array(Inf, dim=c(length(n_h),n_o))
        y.pred<-array(NA,dim=c(nrow(x.test),length(n_h),n_o))
    }else{
        mse.pred<-NULL
        y.pred<-NULL
    }
    # initialize train parameter
    W<-list()
    length(W)<-length(n_h)
    # create randome weights
    if(!is.null(seed)){ # if random seed specified
        set.seed(seed)
    }
    tmp<-runif(max(n_h)*(n_i+1), -1, 1)
    M<-matrix(tmp,nrow=max(n_h),ncol=n_i+1)
    c=1
    for (i in n_h){
        # activation function
        x.with.bias<-cbind(x,rep(1,N))
        o.hidden<-sigmoid(M[1:i,]%*%t(x.with.bias))
        o.hidden.with.bias<-t(rbind(o.hidden,rep(1,N)))

        W[[c]]<-list()
        length(W[[c]])<-n_o
        for (h in 1:n_o){#direct forecast
            W[[c]][[h]]<-ginv(o.hidden.with.bias)%*%y[,h] # based on the first y
            mse.fit[c,h]<-mean((y[,h]-o.hidden.with.bias%*%W[[c]][[h]])^2)
        }

        ### testing
        if(!is.null(x.test)){
            N.t<-nrow(x.test)
            x.test.with.bias<-cbind(x.test,rep(1,N.t))

            o.hidden.t<-sigmoid(M[1:i,]%*%t(x.test.with.bias))
            o.hidden.t.with.bias<-t(rbind(o.hidden.t,rep(1,N.t)))
            for (h in 1:n_o){#direct forecast
                y.pred[,c,h]<-o.hidden.t.with.bias%*%W[[c]][[h]]
                mse.pred[c,h]<-mean((y.test[,h]-y.pred[,c,h])^2)
            }
        }

        c=c+1
    }
    return(list(W=W,M=M,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred))
}

#Pao 1994
# bias on each enhancement node
# direct link from input to output
my.pao.rvfl<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL,seed=NULL){
    y<-as.matrix(y)
    # dim information
    N<-nrow(x)
    n_i<-ncol(x)
    n_o<-ncol(y)
    # initialize test parameter
    mse.fit<-array(Inf, dim=c(length(n_h),n_o))
    if(!is.null(x.test)){
        y.test<-as.matrix(y.test)
        mse.pred<-array(Inf, dim=c(length(n_h),n_o))
        y.pred<-array(NA,dim=c(nrow(x.test),length(n_h),n_o))
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
    c=1 # counter
    for (i in n_h){
        # b, threshold parameter
        a.times.x<-A[1:i,]%*%t(x)
        b<--1*apply(a.times.x,1,mean)
        # enhancement node
        enhancement.node<-t(sigmoid(a.times.x+b))
        # all nodes
        all.node<-cbind(x,enhancement.node)

        W[[c]]<-list()
        length(W[[c]])<-n_o
        for (h in 1:n_o){#direct forecast
            W[[c]][[h]]<-ginv(all.node)%*%y[,h] # based on the first y
            mse.fit[c,h]<-mean((y[,h]-all.node%*%W[[c]][[h]])^2)
        }


        ### testing
        if(!is.null(x.test)){
            a.times.x.t<-A[1:i,]%*%t(x.test)
            # enhancement node
            en.t<-t(sigmoid(a.times.x.t+b))
            # all nodes
            all.node.t<-cbind(x.test,en.t)
            for (h in 1:n_o){#direct forecast
                y.pred[,c,h]<-all.node.t%*%W[[c]][[h]]
                mse.pred[c,h]<-mean((y.test[,h]-y.pred[,c,h])^2)
            }
        }
        c=c+1
    }
    rownames(mse.pred)<-rownames(mse.fit)<-n_h
    colnames(mse.pred)<-colnames(mse.fit)<-1:ncol(y)
    return(list(W=W,A=A,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred))
}

#Chen 1999
# bias on each enhancement node
# direct link from input to output
my.chen.rvfl<-function(x,y,n_h=2*ncol(x),x.test=NULL,y.test=NULL,seed=NULL){
    y<-as.matrix(y)
    # dim information
    N<-nrow(x)
    n_i<-ncol(x)
    n_o<-ncol(y)
    # initialize test parameter
    mse.fit<-array(Inf, dim=c(length(n_h),n_o))
    if(!is.null(x.test)){
        y.test<-as.matrix(y.test)
        mse.pred<-array(Inf, dim=c(length(n_h),n_o))
        y.pred<-array(NA,dim=c(nrow(x.test),length(n_h),n_o))
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
        # b, threshold parameter
        a.times.x<-A[1:i,]%*%t(x)
        b<-runif(i,-1,1)
        # enhancement node
        enhancement.node<-t(tanh(a.times.x+b))
        # all nodes
        all.node<-cbind(x,enhancement.node)

        W[[c]]<-list()
        length(W[[c]])<-n_o
        for (h in 1:n_o){#direct forecast
            W[[c]][[h]]<-ginv(all.node)%*%y[,h] # based on the first y
            mse.fit[c,h]<-mean((y[,h]-all.node%*%W[[c]][[h]])^2)
        }

        ### testing
        if(!is.null(x.test)){
            a.times.x.t<-A[1:i,]%*%t(x.test)
            # enhancement node
            en.t<-t(tanh(a.times.x.t+b))
            # all nodes
            all.node.t<-cbind(x.test,en.t)
            for (h in 1:n_o){#direct forecast
                y.pred[,c,h]<-all.node.t%*%W[[c]][[h]]
                mse.pred[c,h]<-mean((y.test[,h]-y.pred[,c,h])^2)
            }
        }
        c=c+1
    }
    return(list(W=W,A=A,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred))
}

#### cross validation function
cv_fn<-function(x,y,n_h_range,method.text,cv=5){
    # cv to find optimal # of enhancement node
    trn<-list(data=x,labels=y) # put into list
    cv.list<-my_ts_cv_partition(trn,cv=cv)
    cv.mse.grid<-foreach(c=1:cv, .combine='cbind',.export=c('sigmoid','my.pao.rvfl','my.schmidt.rvfl','my.chen.rvfl'),.packages='MASS')%do%{
        x<-cv.list$trn.data[[c]]
        y<-cv.list$trn.labels[[c]]
        x.test<-cv.list$val.data[[c]]
        y.test<-cv.list$val.labels[[c]]

        cv.model<-eval(parse(text=method.text))(x,y,n_h_range,x.test,y.test)
        return(cv.model$mse.pred)
    }
    return(cv.mse.grid)
}

# k-fold CV version of the RVFLs
my.cv.rvfl<-function(x,y,n_h_range,cv=5,methods=c('pao','chen','schmidt'),x.test=NULL,y.test=NULL){
    # determine which method(s) to use
    best.n_h<-method.text<-rep(NA,length(methods))
    cv.mse.grid<-best.model<-list()
    length(cv.mse.grid)<-length(best.model)<-length(methods)
    predict<-array(NA,dim=c(nrow(y.test),length(methods),ncol(y.test)))

    for (i in 1:length(methods)){
        method.text[i]=paste('my.',methods[i],'.rvfl',sep='')
        cv.mse.grid[[i]]<-cv_fn(x,y[,1],n_h_range,method.text[i],cv) # function cv()
        n.cv.mse<-apply(cv.mse.grid[[i]],1,mean)
        best.n_h[i]<-n_h_range[which(n.cv.mse==min(n.cv.mse,na.rm=T))]

        best.model[[i]]<- eval(parse(text=method.text[i]))(x,y,best.n_h[i],x.test,y.test)
        predict[,i,]=best.model[[i]]$predict
    }
    return(list(model=best.model,predict=predict,mse.grid=cv.mse.grid,cv=cv,methods=methods))

}

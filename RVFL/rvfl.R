# define function sigmoid
sigmoid<-function(x){
    x<-as.matrix(x)
    # output is from 0 to 1
    return(1/(1+exp(-x)))
}


# general k-fold CV version of the RVFLs
my.general.cv.rvfl<-function(x,y,n_h_range,cv=5,is.ts=TRUE,x.test=NULL,y.test=NULL,seed=NULL,
                             input.bias=c(TRUE,FALSE),
                             hidden.bias=c(TRUE,FALSE),
                             direct.link=c(TRUE,FALSE),
                             is.scale=FALSE,scale.method='sigma',
                             act.function='sigmoid'){

    # pre allocation
    fit<-array(NA,dim=c(nrow(y),ncol(y),2,2,2))
    best.n_h.grid<-array(NA,dim=c(2,2,2))
    cv.mse.grid<-array(NA,dim=c(length(n_h_range),cv,2,2,2))
    mean.cv.mse<-array(NA,dim=c(length(n_h_range),2,2,2))
    error<-array(NA,dim=c(14,ncol(y),2,2,2))
    # name the dimensions
    dimnames(best.n_h.grid)<-list(c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'))
    dimnames(mean.cv.mse)<-list(n_h_range,c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'))
    dimnames(cv.mse.grid)<-list(n_h_range,1:cv,c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'))
    dimnames(error)<-list(c('MAE','MAPE','MASE','MSE','RMSE','MdAE','MdAPE','RMSPE','RMdSPE','sMAPE','sMdAPE','MRAE','MdRAE','GMRAE'), 1:ncol(y),c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'))

    if(!is.null(x.test)){
        predict<-array(NA,dim=c(nrow(y.test),ncol(y.test),2,2,2))
        dimnames(predict)<-list(1:nrow(y.test),1:ncol(y.test),c('input.T','input.F'), c('hidden.T','hidden.F'),c('direct.T','direct.F'))
    }

    trn.time.start<-proc.time()
    for( i in 1:length(input.bias)){
        has.input.bias=input.bias[i]
        for (j in 1:length(hidden.bias)){
            has.hidden.bias=hidden.bias[j]
            for (k in 1:length(direct.link)){
                has.direct.link=direct.link[k]
                # cv to choose best number of enhancement node
                cv.mse.grid[,,i,j,k]<-general.cv.fn(x,y,n_h_range,cv,is.ts,seed,has.input.bias,has.hidden.bias,has.direct.link,is.scale,scale.method,act.fn)
                if(nrow(cv.mse.grid)==1){#n_h_range is single
                    mean.cv.mse[,i,j,k]<-median(cv.mse.grid[,,i,j,k],na.rm=T)
                }else{
                    mean.cv.mse[,i,j,k]<-apply(cv.mse.grid[,,i,j,k],1,median)
                }
                best.n_h.grid[i,j,k]<-n_h_range[which(mean.cv.mse[,i,j,k]==min(mean.cv.mse[,i,j,k],na.rm=T))]

                # use best.n_h.grid to choose train best model
                best.model<-my.general.rvfl(x,y,best.n_h.grid[i,j,k],x.test,y.test,seed,
                                            has.input.bias,has.hidden.bias,has.direct.link,
                                            is.scale,scale.method,act.fn)
                fit[,,i,j,k]<-best.model$fit[,,1]
                tst.time.start<-proc.time()
                if(!is.null(x.test)){
                    predict[,,i,j,k]<-best.model$predict[,,1]
                    # calculate error measure
                    error[,,i,j,k]<-my_forecasting_measure(predict[,,i,j,k],y.test)
                }
            }
        }
    }
    tst.time<-proc.time()-tst.time.start
    trn.time<-proc.time()-trn.time.start-tst.time

    return(list(error=error,fit=fit,predict=predict,best.n_h=best.n_h.grid,mse.grid=cv.mse.grid,trn.time=trn.time,tst.time=tst.time))

}

#### cross validation function
general.cv.fn<-function(x,y,n_h_range,cv=5,is.ts=TRUE,seed=NULL,
                        has.input.bias=TRUE,
                        has.hidden.bias=TRUE,
                        has.direct.link=TRUE,
                        is.scale=FALSE,scale.method='sigma',
                        act.fn='sigmoid'){
    # cv to find optimal # of enhancement node
    trn<-list(data=x,labels=y) # put into list
    cv.list<-my_cv_partition(trn,cv=cv,is.ts=is.ts)
    mse.grid<-array(NA,dim=c(length(n_h_range),cv))
    for (c in (1:cv.list$k)){ # loop for k fold
        trn.data<-cv.list$trn.data[[c]]
        trn.labels<-cv.list$trn.labels[[c]]
        val.data<-cv.list$val.data[[c]]
        val.labels<-cv.list$val.labels[[c]]
        if(!is.vector(trn.labels)&&!is.factor(trn.labels)){
            trn.labels=trn.labels[,1]
            val.labels=val.labels[,1]
        }
        m<-my.general.rvfl(trn.data,trn.labels,n_h_range,val.data,val.labels,seed,has.input.bias,has.hidden.bias,has.direct.link,is.scale,scale.method,act.fn)
        mse.grid[,c]<-my_forecasting_measure(m$predict[,1,],matrix(rep(val.labels,length(n_h_range)),ncol=length(n_h_range)))['MSE',]

    }

    return(mse.grid)
}


# my general RVFL:
# input bias 0,1
# hidden bias 0,1
# direct link 0,1
# act fn sigm,tanh
my.general.rvfl<-function(x,y,n_h_range=2*ncol(x),x.test=NULL,y.test=NULL,seed=NULL,
                          has.input.bias=TRUE,
                          has.hidden.bias=TRUE,
                          has.direct.link=TRUE,
                          is.scale=TRUE,scale.method='sigma',
                          act.fn='sigmoid',S=1){
    n_h=max(n_h_range)
    y<-as.matrix(y) # convert to row matrix/vector
    x<-as.matrix(x)
    # dim information
    N<-nrow(x) # sample size
    n_i<-ncol(x) #feature size
    n_o<-ncol(y) # horizon
    fit<-array(NA,dim=c(nrow(x),n_o,length(n_h_range)))
    # initialize train parameter
    mse.fit<-array(Inf, dim=c(length(n_h_range),n_o))
    # initialize test parameter
    if(!is.null(x.test)){
        y.test<-as.matrix(y.test)
        x.test<-as.matrix(x.test)
        mse.pred<-array(Inf, dim=c(length(n_h_range),n_o))
        y.pred<-array(NA,dim=c(nrow(x.test),n_o,length(n_h_range)))
    }else{
        mse.pred<-NULL
        y.pred<-NULL
    }

    # create randome weights
    if(!is.null(seed)){ # with pre defined seed
        set.seed(seed)
    }
    tmp<-runif((n_h+1)*(n_i+1), -S, S) # +1 for bias
    #tmp<-rnorm((n_h+1)*(n_i+1), 0, 1) #
    A<-matrix(tmp,nrow=n_h+1,ncol=n_i+1)
    W<-list()
    length(W)<-length(n_h_range)
    counter=1
    for (i in n_h_range){
        print(paste('nh=',i))
        if(has.input.bias){ # has input.bias
            x.new<-cbind(x,rep(S,N))
            RV<-A[1:i,]
        }else{ # no input.bias
            x.new<-x
            RV<-A[1:i,-1] # less one column
        }
        enhancement.node.input<-RV%*%t(x.new)

        if(is.scale){
            if(scale.method=='sigma'){
                scale.factor=2/3
                center.en.input<-mean(enhancement.node.input)
                scale.en.input<-sd(enhancement.node.input)*scale.factor
                scaled.enhancement.node.input<-scale(enhancement.node.input,
                                                     center=rep(center.en.input,N),
                                                     scale=rep(scale.en.input,N))
            }else{ #quantile
                #q<-quantile(enhancement.node.input,c(0.05,0.95))
                q<-apply(enhancement.node.input,1,quantile, c(0.05,0.95))
                s<-c(-log(1/0.05-1),-log(1/0.95-1))
                scaled.enhancement.node.input<-(enhancement.node.input-q[1,])*(s[2]-s[1])/(q[2,]-q[1,])+s[1]
            }
        }else{
            scaled.enhancement.node.input<-enhancement.node.input
        }
        enhancement.node<-eval(parse(text=act.fn))(scaled.enhancement.node.input)
        if(has.hidden.bias){ # has hidden bias
            enhancement.node.new<-t(rbind(enhancement.node,rep(S,N)))
        }else{
            enhancement.node.new<-t(enhancement.node)
        }
        if(has.direct.link){ # direct link
            all.node<-cbind(x.new,enhancement.node.new)
        }else{
            all.node<-enhancement.node.new
        }

        # output layer Least square
        W[[counter]]<-list() # initialize output weights
        length(W[[counter]])<-n_o
        for (h in 1:n_o){#direct forecast
            W[[counter]][[h]]<-ginv(all.node)%*%y[,h] # based on the first y
            fit[,h,counter]<-all.node%*%W[[counter]][[h]]
            mse.fit[counter,h]<-mean((y[,h]-fit[,h,counter])^2)
        }

        ### testing
        if(!is.null(x.test)){
            N.t<-nrow(x.test)
            if(has.input.bias){ # has input.bias
                x.t.new<-cbind(x.test,rep(1,N.t))
                RV.t<-A[1:i,]
            }else{ # no input.bias
                x.t.new<-x.test
                RV.t<-A[1:i,-1] # less one column
            }
            test.enhancement.node.input<-RV.t%*%t(x.t.new)


            if(is.scale){
                if(scale.method=='sigma'){
                    scaled.test.enhancement.node.input<-scale(test.enhancement.node.input,
                                                              center=rep(center.en.input,N.t),
                                                              scale=rep(scale.en.input,N.t))
                }else{ #quantile
                    scaled.test.enhancement.node.input<-(test.enhancement.node.input-q[1,])*(s[2]-s[1])/(q[2,]-q[1,])+s[1]
                }
            }else{
                scaled.test.enhancement.node.input<-test.enhancement.node.input
            }

            en.t<-eval(parse(text=act.fn))(scaled.test.enhancement.node.input)
            if(has.hidden.bias){ # has hidden bias
                en.t.new<-t(rbind(en.t,rep(1,N.t)))
            }else{
                en.t.new<-t(en.t)
            }
            if(has.direct.link){ # direct link
                all.node.t<-cbind(x.t.new,en.t.new)
            }else{
                all.node.t<-en.t.new
            }


            for (h in 1:n_o){#direct forecast
                y.pred[,h,counter]<-all.node.t%*%W[[counter]][[h]]
                mse.pred[counter,h]<-mean((y.test[,h]-y.pred[,h,counter])^2)
            }
        }
        counter=counter+1
    }
    return(list(W=W,A=A,mse.fit=mse.fit,mse.pred=mse.pred,predict=y.pred,fit=fit))
}


# range determination
# input bias: [0 S]
# random weights: [-S S]
# S =2^seq(-5,5,0.5)
# grid search

#### cross validation function
my.general.cv.with.scale.fn<-function(x,y,n_h_range,cv=5,is.ts=TRUE,seed=NULL,
                        has.input.bias=TRUE,
                        has.hidden.bias=FALSE,
                        has.direct.link=TRUE,
                        is.scale=FALSE,scale.method='sigma',
                        S_range=2^seq(-5,5,0.5),
                        act.fn='sigmoid'){
    # cv to find optimal # of enhancement node
    trn<-list(data=x,labels=y) # put into list
    cv.list<-my_cv_partition(trn,cv=cv,is.ts=is.ts)
    mse.grid<-array(NA,dim=c(length(S_range),length(n_h_range),cv))
    for (c in (1:cv.list$k)){ # loop for k fold
        trn.data<-cv.list$trn.data[[c]]
        trn.labels<-cv.list$trn.labels[[c]]
        val.data<-cv.list$val.data[[c]]
        val.labels<-cv.list$val.labels[[c]]
        if(!is.vector(trn.labels)&&!is.factor(trn.labels)){
            trn.labels=trn.labels[,1]
            val.labels=val.labels[,1]
        }
        counter=1
        for (S in S_range){
            print(paste('S=',S))
            m<-my.general.rvfl(trn.data,trn.labels,n_h_range,val.data,val.labels,seed,has.input.bias,has.hidden.bias,has.direct.link,is.scale,scale.method,act.fn, S=S)
            mse.grid[counter,,c]<-my_forecasting_measure(m$predict[,1,],matrix(rep(val.labels,length(n_h_range)),ncol=length(n_h_range)))['MSE',]
            counter=counter+1
        }
    }

    return(mse.grid)
}

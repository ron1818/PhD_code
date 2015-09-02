# my random forests
my_cv_RF<-function(trn,tst,
                is.roll=FALSE,is.reest=FALSE,is.ts=FALSE,classwt=NULL,
                type=if(is.factor(trn$labels)) 'classification' else 'regression'){
  if(type=='classification'){
    is.roll=FALSE # classification cannot roll forecast
  }
  
  if(is.null(ncol(trn$labels))){#only one horizon
    horizon=1
    h.fit<-rep(0,length(trn$labels))
    h.pred<-rep(0,length(tst$labels))
  }else{
    horizon<-ncol(trn$labels)
    h.fit<-array(0,dim=dim(trn$labels))
    h.pred<-array(0,dim=dim(tst$labels))
  }
  
  trn.time.start<-proc.time()
  # rolling forecast
  h.model<-list()
  length(h.model)<-horizon
  
  if(is.roll==TRUE){# rolling forecasing
    # not run for classification
    for (h in 1:horizon){
      # assign labels
      if(horizon==1){ # no need roll, only one data
        h.trn.labels<-trn$labels# 1 step ahead first
        h.tst.labels<-tst$labels# 1 step ahead first
      }else{ # multiple horizon
        h.trn.labels<-trn$labels[,h]# 1 step ahead first
        h.tst.labels<-tst$labels[,h]# 1 step ahead first
      }
      # assign inoput data
      if(h==1){ # first time, always train
        data.width<-ncol(trn$data)
        h.trn.data<-trn$data
        h.tst.data<-tst$data
        # first time, always train
        h.model[[h]]=randomForest(h.trn.data, h.trn.labels, xtest=h.tst.data,ytest=h.tst.labels,
                                  type=type,keep.forest=TRUE,classwt=classwt)$model
      }else{ # h>1, chk reest or not
        # update new dataset for training
        h.trn.data<-cbind(h.trn.data[,-1],h.fit[,h-1])
        h.tst.data<-cbind(h.tst.data[,-1],h.pred[,h-1])
        
        # check if need to retrain randomForest
        if(is.reest){#re estimation of randomForest
          h.model[[h]]=randomForest(h.trn.data,h.trn.labels,xtest=h.tst.data,ytest=h.tst.labels,
                                    type=type,keep.forest=TRUE,classwt=classwt)
        }else{ # copy first model
          h.model[[h]]=h.model[[1]]
        } 
      }# end if h==1
      h.pred[,h]<-predict(h.model[[h]],h.tst.data)
      h.fit[,h]<-predict(h.model[[h]],h.trn.data)
    } # end for horizon
  }else{# not rolling, direct forecasting, bypass is.reest: need to reest by force
    trn.data<-trn$data
    tst.data<-tst$data
    for (h in 1:horizon){
      # assign labels
      if(horizon==1){ # no need roll, only one data
        h.trn.labels<-trn$labels# 1 step ahead first
        h.tst.labels<-tst$labels# 1 step ahead first
      }else{ # multiple horizon
        h.trn.labels<-trn$labels[,h]# 1 step ahead first
        h.tst.labels<-tst$labels[,h]# 1 step ahead first
      }
      # train and test
      h.model[[h]]=randomForest(trn.data, h.trn.labels,xtest=tst.data,ytest=h.tst.labels,
                                type=type,keep.forest=TRUE,classwt=classwt)
      
      if(horizon==1){ # no need roll, only one data
        h.pred<-h.model[[h]]$test$predict
        h.fit<-h.model[[h]]$predict
      }else{ # multiple horizon
        h.pred[,h]<-h.model[[h]]$test$predict
        h.fit[,h]<-h.model[[h]]$predict
      }
    }
    
  }
  trn.time<-proc.time()-trn.time.start
  
  tst.time.start<-proc.time()
  for (h in 1:horizon){
  	# assign labels
  	if(horizon==1){ # no need roll, only one data
  		h.trn.labels<-trn$labels# 1 step ahead first
  		h.tst.labels<-tst$labels# 1 step ahead first
  	}else{ # multiple horizon
  		h.trn.labels<-trn$labels[,h]# 1 step ahead first
  		h.tst.labels<-tst$labels[,h]# 1 step ahead first
  	}
  	# train and test
  	h.model[[h]]=randomForest(trn.data, h.trn.labels,xtest=tst.data,ytest=h.tst.labels,
  														type=type,keep.forest=TRUE,classwt=classwt)
  	
  	if(horizon==1){ # no need roll, only one data
  		h.pred<-h.model[[h]]$test$predict
  		h.fit<-h.model[[h]]$predict
  	}else{ # multiple horizon
  		h.pred[,h]<-h.model[[h]]$test$predict
  		h.fit[,h]<-h.model[[h]]$predict
  	}
  }
  tst.time<-proc.time()-tst.time.start
  
  return(list(model=h.model,predict=h.pred,fit=h.fit,trn.time=trn.time,tst.time=tst.time))
}
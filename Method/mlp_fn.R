# MLP function
require(RSNNS)

my_cv_ann<-function(trn, tst, n_h_range, method='mlp', is.ts=FALSE,
                    type=if(is.factor(trn$labels)) 'classification' else 'regression',
                    cv=5,cv.criteria=NULL,scale=c(0,1),...){
  # scale to [0,1] to fit sigmo
  scaled.trn.data<-my_scale_matrix(trn$data,scale)  
  scaled.tst.data<-my_scale_matrix(tst$data,scale,scaled.trn.data$range)
  if(type=='regression'){
    scaled.trn.labels<-my_scale_matrix(as.matrix(trn$labels),scale)	
    scaled.tst.labels<-my_scale_matrix(as.matrix(tst$labels),scale,scaled.trn.labels$range)
    scaled.trn.labels<-scaled.trn.labels$scaled
    scaled.tst.labels.range<-scaled.tst.labels$range
    scaled.tst.labels<-scaled.tst.labels$scaled
    linOut=TRUE
    if(is.null(cv.criteria))
      cv.criteria='MSE'
  }else{ # classification
    # convert factor to numeric
    scaled.trn.labels<-trn$labels<-decodeClassLabels(trn$labels)
    scaled.tst.labels<-tst$labels<-decodeClassLabels(tst$labels)
    linOut=FALSE
    if(is.null(cv.criteria))
      cv.criteria='F_score'
  }
  
  trn.time.start<-proc.time()
  
  if(cv>1){# cv
    # parallel training to get cv accuracy, # for ts data
    cv.list<-my_cv_partition(list(data=scaled.trn.data$scaled,labels=scaled.trn.labels),
                             is.ts=is.ts,cv=cv)
    measure.grid<-array(0,dim=c(length(n_h_range),cv)) # initialize
    for(i in seq_along(n_h_range)) # loop to choose hidden layer number
    {
      for (j in (1:cv.list$k)){ # loop for k fold
        trn.data<-cv.list$trn.data[[j]]
        trn.labels<-cv.list$trn.labels[[j]]
        val.data<-cv.list$val.data[[j]]
        val.labels<-cv.list$val.labels[[j]]
        
#         if(!is.vector(trn.labels)&&!is.factor(trn.labels)){
#           trn.labels=trn.labels[,1]
#           val.labels=val.labels[,1]
#         }
        m<-eval(parse(text=method))(trn.data, trn.labels, 
                                    inputsTest=val.data, targetsTest=val.labels,
                                    size=c(n_h_range[i]),linOut=linOut)
        pred<-m$fittedTestValues
        if(type!='classification'){
          measure.grid[i,j]<-my_forecasting_measure(pred,val.labels)[cv.criteria,1]
        }else{
          pred<-encodeClassLabels(pred)
          val.labels<-encodeClassLabels(val.labels)
          # reciprocal, small is better
          measure.grid[i,j]<-1/my_classification_measure(pred,val.labels)[[cv.criteria]]
        }
      }
    }
    
    # average out to find cv accuracy
    measure.score<-rowMeans(measure.grid,na.rm=T)
    # find best param from the grid
    best.param.idx<-which(measure.score==min(measure.score),arr.ind=T)
    best.n_h<-n_h_range[best.param.idx]
  }else{ #no cv
    best.n_h<-n_h_range[1]
    best.param.idx<-1
    measure.grid<-NULL
  }
  
  # test
  trn.time<-proc.time()-trn.time.start
  
  tst.time.start<-proc.time()
  model<-eval(parse(text=method))(scaled.trn.data$scaled, scaled.trn.labels, 
             inputsTest=scaled.tst.data$scaled, targetsTest=scaled.tst.labels,
             size=c(best.n_h),linOut=linOut)
  
  # scale back
  predict=model$fittedTestValues
  if(type=='regression'){
    predict<-my_scale_matrix(predict,scaled.tst.labels.range)$scaled
    tst.labels<-tst$labels
  }else{# classification
    predict<-encodeClassLabels(predict)
    tst.labels<-encodeClassLabels(scaled.tst.labels)
  
  }
  tst.time<-proc.time()-tst.time.start
  return(list(model=model,predict=predict,tst.labels=tst.labels,n_h=best.n_h, measure.grid=measure.grid, trn.time=trn.time,tst.time=tst.time))
  
}
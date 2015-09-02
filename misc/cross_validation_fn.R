# n-fold cv
my_cv_partition<-function(trn,is.ts=FALSE,cv=5){
  trn.length=nrow(trn$data)
  trn.width=ncol(trn$data)
  if(is.ts){ # for ts, sequential
    partition.number=cv*2-1 # rolling cv
    # make factor to partition
    f<-rep(1:partition.number,each=floor(trn.length/partition.number))
  }else{ # non cv partition
  # make factor to partition
    f<-ceiling(runif(trn.length,0,cv)) # cv=1 means no cv
  }
  # initialize 
  x<-y<-list()
  x.labels<-y.labels<-list()
  length(x)<-length(y)<-length(x.labels)<-length(y.labels)<-cv

  if(cv>1){# valid cv number
    for (c in 1:cv){
      y[[c]]<-trn$data[which(f==c),] # validation data
      x[[c]]<-trn$data[which(f!=c),] # training data
      if(is.null(nrow(trn$labels))){
        y.labels[[c]]<-trn$labels[which(f==c)] # validataion labels
        x.labels[[c]]<-trn$labels[which(f!=c)] # training labels
      }else{
        y.labels[[c]]<-trn$labels[which(f==c),,drop=F] # validataion labels
        x.labels[[c]]<-trn$labels[which(f!=c),,drop=F] # training labels
      }
    }
    return(list(trn.data=x, trn.labels=x.labels, val.data=y, val.labels=y.labels,k=cv))
  }else{#no cv
    return(NULL)
  }
}

# n-fold cv for time series
# rolling
my_ts_cv_partition<-function(trn,cv=5){
  #trn<-list(data,labels)
  trn.length=nrow(trn$data)
  trn.width=ncol(trn$data)
  partition.number=cv*2-1 # rolling cv
  # make factor to partition
  f<-rep(1:partition.number,each=floor(trn.length/partition.number))

  # initialize 
  x<-y<-list()
  x.labels<-y.labels<-list()
  length(x)<-length(y)<-length(x.labels)<-length(y.labels)<-cv
  
  if(cv>1){# valid cv number
    for (c in (1:cv)){
      x[[c]]<-trn$data[which(f>=c & f<=c+(cv-2)),]# training data
      x.labels[[c]]<-as.matrix(trn$labels)[which(f>=c & f<=c+(cv-2)),]# training labels
      y[[c]]<-trn$data[which(f==cv+c-1),]# validation data
      y.labels[[c]]<-as.matrix(trn$labels)[which(f==cv+c-1),]# validataion labels
    }
    return(list(trn.data=x, trn.labels=x.labels, val.data=y, val.labels=y.labels,k=cv))
  }else{#no cv
    return(NULL)
  }
}

# 5X2 CV F test
my_52cv<-function(acc1, acc2, alpha){
  n_folds=5
  n_cv=2
  # initialize
  p<-array(0, dim=c(n_folds, n_cv))
  p.bar<-array(0, dim=c(n_folds, 1))
  s.sqr<-array(0, dim=c(n_folds, 1))
  
  # calculate p and s.sqr
  for (i in (1:n_folds))
  {
    for (j in (1:n_cv))
    {
      p[i,j]<-acc1[i,j]-acc2[i,j]
    }
    p.bar[i]<-0.5*(p[i,1]+p[i,2])
    s.sqr<-(p[i,1]-p.bar[i])^2+(p[i,2]-p.bar[i])^2
  }
  # calculate f
  f<-sum(p^2)/(2*sum(s.sqr))
  f.critial<-qf(1-alpha, 10, 5)
  
  # The critical value is the number that the test statistic must exceed to reject the test
  
  if (f > f.critial)
  {
    if(mean(acc1)>mean(acc2))
    {
      return(list(decision=1, f.value=f)) # acc1 sig. better than acc2
    }else
    {
      return(list(decision=0, f.value=f)) # acc1 acc2 no evidence differ
    }
  }else
  {
    return(list(decision=0, f.value=f)) # acc1 acc2 no evidence differ
  }
}
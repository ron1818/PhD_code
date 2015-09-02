# convert ts into vectors with LAG as the col.width
my_ts2vector<-function(ts,LAG)
{
  # training data
  x.width<-LAG
  x.length<-length(ts)-x.width+1
  x<-array(0,dim=c(x.length, x.width))
  
  for (i in (1:x.width))
  {
    x[,i]=ts[i:(i+x.length-1)]
  }
  return(x)
}

# scale ts to within range
my_scale_ts<-function(ts, scaled.range=c(0,1),ts.range=range(ts))
{
  ts.max<-max(ts.range)
  ts.min<-min(ts.range)
  d.max<-max(scaled.range)
  d.min<-min(scaled.range)
  
  scaled.ts<-(ts-ts.min)*(d.max-d.min)/(ts.max-ts.min)+d.min
  return(list(scaled=scaled.ts, range=c(ts.min,ts.max)))
}


# scale ts to within range
my_scale_matrix<-function(mat, scaled.range=c(-1,1),mat.range=NULL)
{
	if(is.null(dim(scaled.range))){# single range for all column
		scaled.range=matrix(scaled.range,2,ncol(mat))	 # repeat
	}
	# find original max and min
	if(is.null(mat.range)){
		mat.max<-apply(mat,2,max,na.rm=T)
		mat.min<-apply(mat,2,min,na.rm=T)
	}else{
		mat.max=mat.range[2,]
		mat.min=mat.range[1,]
	}
	
	scaled.mat<-matrix(NA,nrow(mat),ncol(mat))
	for(i in 1:ncol(mat)){
		vector<-mat[,i] # column by column
		s.range<-scaled.range[,i]
		vector.max<-max(vector, na.rm=T)
		vector.min<-min(vector, na.rm=T)
		d.max<-max(s.range)
		d.min<-min(s.range)
		scaled.mat[,i]<-(vector-vector.min)*(d.max-d.min)/(vector.max-vector.min)+d.min
	}
	return(list(scaled=scaled.mat,range=rbind(mat.max,mat.min)))
}

# split vector with ratio (training data/total data)
my_vector_split<-function(ts,HORIZON=1,LAG=2,RATIO=0.7,is.point=F)
{
	if(is.point){#point forecasting, need to delect intermediate data
		org.HORIZON=HORIZON
		HORIZON=max(HORIZON)
	}
  # training data
  x.width<-LAG+HORIZON
  x.length<-length(ts)-x.width+1
  x<-array(0,dim=c(x.length, x.width))
  
  for (i in (1:x.width))
  {
    x[,i]=ts[i:(i+x.length-1)]
  }
  
  historical<-x[,(1:LAG),drop=F]
  future<-x[,LAG+(1:HORIZON),drop=F]
  
  # split according to ratio
  # trn.idx<-1:round(nrow(historical)*RATIO)
	tst.length<-length(ts)-round(length(ts)*RATIO)
	tst.idx<-(nrow(historical)-tst.length+1):nrow(historical)
  
  trn.data<-as.matrix(historical)[-tst.idx,,drop=F]
  trn.labels<-as.matrix(future)[-tst.idx,,drop=F]
	if(tst.length==0){# no tst
		tst.data=NULL
		tst.labels=NULL
	}else{
		tst.data<-as.matrix(historical)[tst.idx,,drop=F]
		tst.labels<-as.matrix(future)[tst.idx,,drop=F]
	}
	
	if(is.point){ # discard unused points
		trn.labels<-trn.labels[,org.HORIZON,drop=F]
		if(tst.length!=0)
			tst.labels<-tst.labels[,org.HORIZON,drop=F]
	}
  
  trn<-list(data=trn.data,labels=trn.labels)
  tst<-list(data=tst.data,labels=tst.labels)
  
  return(list(trn=trn,tst=tst))
}
# calculate distance
my_knn_dist<-function(x,y,w=rep(TRUE,ncol(y)),k=round(nrow(x)/2)){
	# distance matrix of y vs. x
	y.dist<-array(0, dim=c(nrow(y), nrow(x)))
	knn.idx<-knn.dist<-array(0, dim=c(nrow(y), k))
	for (i in (1:nrow(y)))
	{
		# duplicate ith row of y
		ith.y<-array(rep(as.numeric(y[i,]),each=nrow(x)),dim=dim(x))
		# weighted euclidean distance
		y.dist[i,]<-sqrt((ith.y-x)^2 %*% w)
		# select first nearest k indices
		knn.idx[i,]<-order(y.dist[i,], decreasing=F)[1:k]
		knn.dist[i,]<-y.dist[i,knn.idx[i,]]
	}
	return(list(dist.mat=y.dist,knn.idx=knn.idx,knn.dist=knn.dist))
}

# rmse k vs. f for a given set
my_knn_rmse_grid<-function(trn.data,trn.labels,tst.data,tst.labels,k_range,f_range,w){
	y<-tst.data # testing data
	y.labels<-tst.labels # testing labels
	x<-trn.data # training data
	x.labels<-trn.labels # training labels
	# initialize
	rmse.grid<-array(Inf,dim=c(nrow(f_range),length(k_range)))
	# loop on f
	for (i in 1:nrow(f_range)){
		f=f_range[i,]
		fth.knn.model<-my_knn_dist(x,y,f,max(k_range)) # find model
		fth.knn.idx<-fth.knn.model$knn.idx
		fth.knn.dist<-fth.knn.model$knn.dist
		for (j in seq_along(k_range)){
			k=k_range[j]
			kth.idx<-fth.knn.idx[,1:k]
			kth.dist<-fth.knn.dist[,1:k]
			knn.pred<-array(x.labels[kth.idx],dim=c(dim(y)[1], k))
			# linear weight decay
			norm.weight<-(1/kth.dist)/apply((1/kth.dist),1,sum)
			y.pred<-apply(knn.pred*norm.weight,1,sum)
			rmse.grid[i,j]=sqrt(sum((y.labels-y.pred)^2))
		}
	}
	return(rmse.grid)
}
##########END Functions###########


my_knn<-function(trn,tst,k_range=c(3,9,21),f_range=matrix(rep(TRUE,ncol(trn$data)),nrow=1),
								 feature.weight=rep(1,ncol(trn$data)),cv=5,is.ts=FALSE,plot.rmse=FALSE)
{
  # separate data and labels
  x<-trn$data
  x.labels<-trn$labels
  y<-tst$data
  y.labels<-tst$labels
  
  # calculate feature weight
  feature.w<-feature.weight
  # check if need cv or not
  if (cv>1){ # require cv
  	
  	# parallel training to get cv accuracy, # is.ts=TRUE for ts data
  	cv.list<-my_cv_partition(trn,is.ts=is.ts,cv=cv)
  	
  	rmse.grid<-array(NA,dim=c(nrow(f_range),length(k_range),cv))
  	cv.rmse.grid<-array(0,dim=c(nrow(f_range),length(k_range)))
  		for (j in (1:cv.list$k)){ # loop for k fold
  			trn.data<-cv.list$trn.data[[j]]
  			trn.labels<-cv.list$trn.labels[[j]]
  			val.data<-cv.list$val.data[[j]]
  			val.labels<-cv.list$val.labels[[j]]
  			if(!is.vector(trn.labels)&&!is.factor(trn.labels)){
  				trn.labels=trn.labels[,1]
  				val.labels=val.labels[,1]
  			}
  			rmse.grid[,,j]<-my_knn_rmse_grid(trn.data,trn.labels,val.data,val.labels,k_range=k_range,f_range=f_range,feature.w)
  			cv.rmse.grid<-cv.rmse.grid+rmse.grid[,,j]
  		}	
  	# average out to find cv accuracy
  	cv.rmse.grid<-cv.rmse.grid/cv

  	# plot RMSE grid
    if(plot.rmse==TRUE){ # need to plot rmse surface
      z<-cv.rmse.grid[f_range,k_range]
      nrz <- length(f_range)
      ncz <- length(k_range)
      # surface plot
      jet.colors <- colorRampPalette( c("blue", "green") )
      # Generate the desired number of colors from this palette
      nbcol <- 100
      color <- jet.colors(nbcol)
      # Compute the z-value at the facet centres
      zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
      # Recode facet z-values into color indices
      facetcol <- cut(zfacet, nbcol)
      
      persp(f_range,k_range,z, theta=60,phi=30,xlab='f',ylab='k',zlab='rmse',col=color[facetcol],ticktype='detailed',main='RMSE k v.s. f')
    }
  	
  	# find best parameters
    best.param<-which(cv.rmse.grid==min(cv.rmse.grid,na.rm=T),arr.ind=T) 
    f=f_range[best.param[1],]
    k=k_range[best.param[2]]
    }
  else{# no cv, k_range==1, just calculate
    k=k_range[1]
    f=rep(TRUE, ncol(x))
  }
  
  knn.model<-list() # initialize knn model
  y.pred<-array(0,dim=dim(y.labels))
    for ( j in 1:ncol(y.labels)){ # different horizon
      knn.model<-my_knn_dist(x,y,w=f*feature.w,k=k)
      knn.idx<-knn.model$knn.idx
      knn.dist<-knn.model$knn.dist
      norm.weight<-(1/knn.dist)/apply((1/knn.dist),1,sum)
      knn.pred<-array(x.labels[knn.idx,j],dim=c(nrow(y), k))
      y.pred[,j]<-apply(knn.pred*norm.weight,1,sum)
    }
  
  return(list(predict=y.pred,model=knn.model,feature.weight=feature.weight,best.param=best.param))
}


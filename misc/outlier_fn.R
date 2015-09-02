# outlier detection and removal functions
#http://stats.stackexchange.com/questions/1142/simple-algorithm-for-online-outlier-detection-of-a-generic-time-series
my_res_iqr_outliers <- function(x,iqr.factor=1.5,plot=TRUE)
{
  x <- as.ts(x)
  if(frequency(x)>1){ # may have seasonal
    resid <- stl(x,s.window="periodic",robust=TRUE,na.action='na.exclude')$time.series[,3]
  }else{
    tt <- 1:length(x)
    resid <- residuals(loess(x ~ tt, na.action='na.exclude')) # using loess to fit
  }
  resid.q <- quantile(resid,prob=c(0.25,0.75), na.rm=T) # get IQR
  iqr <- diff(resid.q)
  limits <- resid.q + iqr.factor*iqr*c(-1,1) # outside 3IQR => outliers
  score <- abs(pmin((resid-limits[1])/iqr,0) + pmax((resid - limits[2])/iqr,0)) # +ve for outlier, otherwise non outlier
  outlier.idx<-which(score>0)
  if(plot)
  {
    plot(x)
    x2 <- ts(rep(NA,length(x)))
    x2[outlier.idx]<-x[outlier.idx]
    tsp(x2) <- tsp(x)
    points(x2,pch=19,col="red")
    return(invisible(outlier.idx))
  }
  else
    return(outlier.idx)
}

my_window_mad_outliers<-function(x, window=30, threshold=5)
{
  x <- as.ts(x)
  ut<-function(x,threshold){
    median(x, na.rm=T)+threshold*mad(x, na.rm=T)
  }
  z<-rollapply(x, window, ut, threshold=threshold, align='right')
  z <- c(rep(z[1], window-1), z) # Use z[1] throughout the initial period
  outlier <- x > z
  outlier.idx<-which(outlier==T)  
  plot(x, type="l")
  x2 <- ts(rep(NA,length(x)))
  x2[outlier.idx]<-x[outlier.idx]
  tsp(x2) <- tsp(x)  
  points(x2, pch=19,col='red')
  return(outlier.idx)
}

### rearrange data to vector
my_delete_outlier<-function(x,outlier.idx=NULL)
{
  x<-na.approx(x)
  
  if(is.xts(x)||is.ts(x)) timestamp<-time(x) else timestamp=NULL
  
  x.nooutlier<-as.numeric(x)
  
  # using SMA to remove outlier, causal
  for (i in outlier.idx)
  {
    if(i>2)
      x.nooutlier[i]<-0.5*(x.nooutlier[i-2]+x.nooutlier[i-1])
    else if(i<=2&&i>1)
      x.nooutlier[i]<-data(x.nooutlier[i-1])
    else #i==1
      x.nooutlier[i]<-0
  }
  if(!is.null(timestamp))
    x.nooutlier<-xts(x.nooutlier,timestamp)
  
  return(x.nooutlier)
}

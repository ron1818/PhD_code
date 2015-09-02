my_majority_vote<-function(x,class=unique(as.numeric(x)),levels=unique(as.numeric(x)),labels=levels,weight=NULL,byrow=TRUE){
	if(byrow==FALSE) # column of x is sample, row of x is votes
		x=t(x) # need to transpose
	
	if(is.null(weight)) #uniform weighted vote
		weight=rep(1,ncol(x))
	
	weight=rep(weight,nrow(x))
	weight=matrix(weight,nrow(x),byrow=TRUE)
	
	# find votes
	votes<-array(NA,dim=c(nrow(x),length(class)))
	for (i in seq_along(class)){# for different class, calculate weighted votes
		tmp.weight=weight
		tmp.weight[which(x!=class[i])]=0
		votes[,i]=apply(tmp.weight,1,sum)
	}
	colnames(votes)<-class
	
	# find final decision
	agg.x<-rep(NA,nrow(x))
	for (idx in 1:nrow(x)){
		highest.votes<-which(votes[idx,]==max(votes[idx,]))
		agg.x[idx]<-colnames(votes)[highest.votes]
	}

		agg.x=factor(agg.x,labels=labels,levels=levels)
	
	return(list(votes=votes,value=agg.x))
}





### Error Measure
# R. J. Hyndman and A. B. Koehler 2005
my_forecasting_measure<-function(y.hat, y, err.benchmk=y[,1])
{
  y<-as.matrix(y)
  y.hat<-as.matrix(y.hat)
  err.benchmk<-as.matrix(err.benchmk)
  n<-nrow(y)
  
  #Scale dependent: MAE RMSE MAE MdAE
  error<-y-y.hat
  error<-as.matrix(error)
  MSE<-apply(error^2,2,mean, na.rm=T)
  RMSE<-sqrt(MSE)
  MAE<-apply(abs(error),2,mean, na.rm=T)
  MdAE<-apply(abs(error),2,median, na.rm=T)
  
  #Percentage errors: MAPE MdAPE RMSPE RMdSPE
  per.err<-error/y*100
  MAPE<-apply(abs(per.err),2,mean, na.rm=T)
  MdAPE<-apply(abs(per.err),2,median, na.rm=T)
  RMSPE<-sqrt(apply(per.err^2,2,mean, na.rm=T))
  RMdSPE<-sqrt(apply(per.err^2,2,median, na.rm=T))
  
  #Symmetric
  s.per.err<-200*abs(y-y.hat)/(y+y.hat)
  s.per.err<-as.matrix(s.per.err)
  sMAPE<-apply(s.per.err,2,mean, na.rm=T)
  sMdAPE<-apply(s.per.err,2,median, na.rm=T)
  
  #Relative: MRAE MdRAE GMRAE
  rela.err<-error/rep(err.benchmk,ncol(error))
  MRAE<-apply(abs(rela.err),2,mean, na.rm=T)
  MdRAE<-apply(abs(rela.err),2,median, na.rm=T)
  GMRAE<- exp(apply(log(abs(rela.err)),2,mean, na.rm=T))
  
  #Scaled
  y<-as.matrix(y)
  n<-nrow(y)
  MASE<-apply(abs(error),2,sum,na.rm=T)/(n/(n-1)*apply(abs(diff(y)),2,sum,na.rm=T))
  
  error.measure<-rbind(MAE, MAPE, MASE,
               MSE, RMSE,MdAE,
               MdAPE, RMSPE, RMdSPE,
               sMAPE, sMdAPE,
               MRAE, MdRAE, GMRAE)
  rownames(error.measure)<-c('MAE', 'MAPE', 'MASE',
                             'MSE', 'RMSE','MdAE',
                             'MdAPE', 'RMSPE', 'RMdSPE',
                             'sMAPE', 'sMdAPE',
                             'MRAE', 'MdRAE', 'GMRAE')
  colnames(error.measure)<-1:ncol(y)
  return(error.measure)
}

my_forecasting_measure_solar<-function(y.hat, y, err.benchmk=y[,1],sunrise)
{
	sunrise[sunrise==0]=NA
	y<-as.matrix(y)
	y.hat<-as.matrix(y.hat)
	err.benchmk<-as.matrix(err.benchmk)
	n<-nrow(y)
	
	#Scale dependent: MAE RMSE MAE MdAE
	error<-y-y.hat
	error<-as.matrix(error)
	nanerror<-as.matrix(error*sunrise)
	MSE<-apply(nanerror^2,2,mean, na.rm=T)
	RMSE<-sqrt(MSE)
	MAE<-apply(abs(nanerror),2,mean, na.rm=T)
	MdAE<-apply(abs(nanerror),2,median, na.rm=T)
	
	#Percentage errors: MAPE MdAPE RMSPE RMdSPE
	per.err<-error/y*100
	nanper.err<-per.err*sunrise
	MAPE<-apply(abs(nanper.err),2,mean, na.rm=T)
	MdAPE<-apply(abs(nanper.err),2,median, na.rm=T)
	RMSPE<-sqrt(apply(nanper.err^2,2,mean, na.rm=T))
	RMdSPE<-sqrt(apply(nanper.err^2,2,median, na.rm=T))
	
	#Symmetric
	s.per.err<-200*abs(y-y.hat)/(y+y.hat)
	s.per.err<-as.matrix(s.per.err)
	sMAPE<-apply(s.per.err,2,mean, na.rm=T)
	sMdAPE<-apply(s.per.err,2,median, na.rm=T)
	
	#Relative: MRAE MdRAE GMRAE
	rela.err<-error/rep(err.benchmk,ncol(error))
	MRAE<-apply(abs(rela.err),2,mean, na.rm=T)
	MdRAE<-apply(abs(rela.err),2,median, na.rm=T)
	GMRAE<- exp(apply(log(abs(rela.err)),2,mean, na.rm=T))
	
	#Scaled
	y<-as.matrix(y)
	n<-nrow(y)
	MASE<-apply(abs(error),2,sum,na.rm=T)/(n/(n-1)*apply(abs(diff(y)),2,sum,na.rm=T))
	
	error.measure<-rbind(MAE, MAPE, MASE,
											 MSE, RMSE,MdAE,
											 MdAPE, RMSPE, RMdSPE,
											 sMAPE, sMdAPE,
											 MRAE, MdRAE, GMRAE)
	rownames(error.measure)<-c('MAE', 'MAPE', 'MASE',
														 'MSE', 'RMSE','MdAE',
														 'MdAPE', 'RMSPE', 'RMdSPE',
														 'sMAPE', 'sMdAPE',
														 'MRAE', 'MdRAE', 'GMRAE')
	colnames(error.measure)<-1:ncol(y)
	return(error.measure)
}

my_classification_measure<-function(predict,target){
	# 20) error metrics
	# 20.1) The kappa coefficient (Cohen 1960) is a statistical measure of inter-rater agreement for qualitative (categorical) events
	# k=(pa-pe)/(1-pe)
	
	# 20.2) Contingency table
#   contingency_table<-array(NA,dim=c(2,2))
#   rownames(contingency_table)<-colnames(contingency_table)<-c(-1,1)
	tmp<-table(predict, target)
#   if(dim(tmp)==dim(contingency_table))
    contingency_table<-tmp
#   else
    
	TP<-contingency_table[2,2]
	TN<-contingency_table[1,1]
	FP<-contingency_table[2,1]
	FN<-contingency_table[1,2]
	
	ACC<-(TP+TN)/(TP+TN+FP+FN)
	Precision<-TP/(TP+FP)
	Recall<-Sensitivity<-H_rate<-TP/(TP+FN) # also known as probability of detection, sensitivity, hit rate, true pos rate
	Specificity<-TN/(TN+FP) # true nega rate
	F_score<-2*(Precision*Recall)/(Precision+Recall)
	MCC<-(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
	CSI<-TP/(TP+FN+FP) # critical success index (CSI)
	Bias_score<-(TP+FP)/(TP+FN)
	F_rate<-FP/(FP+TN) # false alarm rate
	SR<-1-F_rate # success ration
	KSS<-H_rate-F_rate # Hanssen & Kuipers skill score
	if(is.null(nrow(target)))
		n<-length(target)
	else
		n<-nrow(target)
	
	EDS<-(2*log(TP+FN)/n)/log(TP/n)-1# extreme dependency score
	OR<- (H_rate/(1-H_rate))/(F_rate/(1-F_rate)) # odds ratio
	ORSS<- (OR-1)/(OR+1) # odds ratio skill score
	return(list(contigency_table=contingency_table,Precision=Precision,Recall=Recall,ACC=ACC,Specificity=Specificity,
							F_score=F_score,CSI=CSI,KSS=KSS,EDS=EDS,ORSS=ORSS,SR=SR))
}



my_Friedman<-function(x, alpha=0.05)
{
  N<-dim(x)[1]
  k<-dim(x)[2]
  rank<-x 
  # calculate average rank
  for (i in 1:N)
  {
    rank[i,]<-rank(x[i,], ties.method='average')
  }
  ave_rank<-colMeans(rank)
  # calculate chi square statistics
  chi2<-12*N/(k*(k+1))*(sum(ave_rank^2)-k*(k+1)^2/4)
  # calculate F statistic
  Ff<-(N-1)*chi2/(N*(k-1)-chi2)
  # F statistic theoritical
  dof<-c((k-1),(k-1)*(N-1))
  F_critical<-qf(1-alpha, dof[1], dof[2])
  p.value<-pchisq(chi2,dof[1],lower.tail=FALSE)
  # is_null_reject
  return(list(n0.rej=Ff>F_critical, ave_rank=ave_rank, dof=dof,chi2=chi2,Ff=Ff,p.value=p.value))
}

my_Nemenyi<-function(x,alpha=0.05)
{
	N<-dim(x)[1]
	k<-dim(x)[2]
	#q.alpha
	q.alpha<-qtukey(1-alpha,k,Inf)/sqrt(2)
	# critical difference
	CD<-q.alpha*sqrt(k*(k-1)/(6*N))
	rank<-x 
	# calculate average rank
	for (i in 1:N)
	{
		rank[i,]<-rank(x[i,], ties.method='average')
	}
	ave_rank<-colMeans(rank)
	ave_rank<-sort(ave_rank)
	
	dist.mat<-abs(matrix(rep(ave_rank,k),nrow=k)
                -matrix(rep(ave_rank,k),nrow=k,byrow=TRUE))
  dist.mat.tf<-upper.tri(dist.mat)&(dist.mat<=CD)
  dist.pair<-which(dist.mat.tf,arr.ind=TRUE)
  stack<-array(NA,dim(dist.pair))
	# unique row.idx
  counter=1
  row.uniq<-unique(dist.pair[,1])
  for (i in row.uniq){
    stack[counter,]<-c(i,max(dist.pair[dist.pair[,1]==i,2]))
    counter=counter+1
  }
  dist.pair.row.uniq<-stack[1:(counter-1),,drop=FALSE]
  # unique col.idx
  stack<-array(NA,dim(dist.pair.row.uniq))
  counter=1
  col.uniq<-unique(dist.pair.row.uniq[,2])
	for (i in col.uniq){
	  stack[counter,]<-c(min(dist.pair.row.uniq[dist.pair.row.uniq[,2]==i,1]),i)
	  counter=counter+1
	}
  dist.pair.final<-stack[1:(counter-1),,drop=FALSE]
  
	# visualization
	require(ggplot2)
	require(grid)
	df.ave_rank<-as.data.frame(ave_rank)
	p<-ggplot(df.ave_rank,aes(x=ave_rank,y=0))+geom_point(size=4)+geom_text(label=rownames(df.ave_rank), size=4,vjust=0,hjust=1.5,angle=90)
  
	p1<-p+annotate('segment',size=1,x=ave_rank[dist.pair.final[,1]],xend=ave_rank[dist.pair.final[,2]],y=seq(0.1,0.2,length.out=nrow(dist.pair.final)),yend=seq(0.1,0.2,length.out=nrow(dist.pair.final)))+ylim(-1,0.3)+annotate('segment',x=1,xend=1+CD,y=-0.8,yend=-0.8,arrow=arrow(ends="both", angle=90,length=unit(6,'pt')))+annotate('text',x=1+CD/2,y=-0.8+0.1,label=paste('CD',round(CD,3),sep='='))+ylab('')+xlab('Average Rank')+theme_bw()+theme(axis.text.y = element_blank(),axis.ticks=element_blank())

  return(list(p=p1,ave_rank=ave_rank,dist.pair=dist.pair.final, CD=CD))

}

my_Wilcoxon<-function(x, y)
{
  N<-length(x)
  diff<-x-y
  rank<-rank(abs(diff), ties.method='average')
  sign.rank<-sign(diff)
  
  r.plus<-sum(rank*as.numeric(sign.rank>=0))
  r.minus<-sum(rank*as.numeric(sign.rank<=0))
  T<-min(r.plus, r.minus)
  z<-(T-0.25*N*(N+1))/sqrt(1/24*N*(N+1)*(2*N+1))
  # is_null_reject
  return(list(n0.rej=ifelse(abs(z)>1.96,TRUE,FALSE), z=z))
}

my_Residual_Analysis<-function(x,Lag=15)
{
  if(is.numeric(x))
    N<-length(x)
  else
    N<-nrow(x)
  
  df<-data.frame(idx=1:N,res=x)
  require(ggplot2)
  # plot residuals
  p1<-ggplot(df,aes(x=idx,y=res))+geom_line()+xlab('Time')+ylab('Residuals')
  
  # plot acf
  conf.level <- 0.95
  ciline <- qnorm((1 - conf.level)/2)/sqrt(N)
  bacf <- acf(x, lag.max=Lag,plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  p2 <- ggplot(data = bacfdf, aes(x = lag, y = acf))+geom_point()+
    geom_hline(aes(yintercept =ciline),linetype='dashed',colour='blue') +
    geom_hline(aes(yintercept =-ciline),linetype='dashed',colour='blue') +
    geom_segment(aes(xend = lag, yend = 0))+xlab('Lag')+ylab('ACF')
  
  # plot spectrum
  spectrum.x<-spectrum(x,plot=FALSE)
  spectrum.x.df<-data.frame(freq<-spectrum.x$freq,spec<-spectrum.x$spec)
  p3<-ggplot(data=spectrum.x.df,aes(x=freq,y=spec))+geom_line()+scale_y_log10()+xlab('Frequency')+ylab('Spectrum')
  
  # plot Ljung-box
  p<-rep(NA,Lag)
  for (i in 1:Lag){
    p[i]<-Box.test(x, i, type = "Ljung-Box")$p.value
  }
  p.df<-data.frame(p=p,lag=1:Lag)
  
  p4<-ggplot(data=p.df,aes(x=lag,y=p))+geom_point(size=4,shape=21)+geom_hline(aes(yintercept =0.05),linetype='dashed',colour='blue')+xlab('Lag')+ylab('p-value')+ylim(0,1)
  
  grid.arrange(p1,p2,p3,p4,nrow=4)
}
#### how to cite  #####
# @INPROCEEDINGS{Ren2013CIEL,
#                author = {Ye Ren and Le Zhang and P. N. Suganthan},
#                title = {K{NN} based bagging {SVM} pruning},
#                booktitle = {Proc. IEEE Symposium on Computational Intelligence and Ensemble Learning
#                             (CIEL'13)},
#   year = {2013},
#   address = {Singapore},
#   month = apr
# }



### 0) initialize ###
# setwd('./CIEL2013') # current working 
# required packages, if not installed, please use:
# install.packages() function to install the packages
require(xts)
require(e1071)
require(bootstrap)
require(RSNNS)
require(timeDate)
require(forecast)
# source user defined functions
source('../misc/pre_processing_fn.R')
source('../misc/post_processing_fn.R')
source('../misc/cross_validation_fn.R')
source('../Method/svm_fn.R')
source('../Method/knn_fn.R')

### 1) load data ###
breast.wisconsin<-read.csv('breast_wisconsin.csv')
labels<-breast.wisconsin[,'class'] #2 bengin, 4 malignant
labels[which(labels==2)]=-1 #bengin
labels[which(labels==4)]=1 #maligant
breast.wisconsin[,'class']<-labels # convert to -1, +1

# find complete cases
breast.wisconsin<-breast.wisconsin[complete.cases(breast.wisconsin),]
# scale data
data<-scale(breast.wisconsin[,2:10]) # avoid id
labels<-breast.wisconsin[,'class']

# split to training and testing
RATIO=0.7
trn.number<-round(nrow(data)*RATIO)
trn.idx<-sample(nrow(data),trn.number)
tst.idx<--1*trn.idx

trn.data<-data[trn.idx,]
trn.labels<-labels[trn.idx]
tst.data<-data[tst.idx,]
tst.labels<-labels[tst.idx]

# bootstrap
B=100
idx.length<-nrow(trn.data)
boot<-bootstrap(1:idx.length, B, as.numeric)
idx.boot<-boot$thetastar # oob.idx=-idx.boot

# training
model<-list()
prediction<-probability<-array(NA,c(nrow(trn.data),B))
for ( i in 1:B){
	model[[i]]<-svm(trn.data[idx.boot[,i],],trn.labels[idx.boot[,i]],type='C-classification',probability=TRUE)
	tmp<-predict(model[[i]],trn.data,probability=TRUE)
	prediction[,i]<-as.numeric(levels(tmp)[tmp])
	probability[,i]<-apply(attr(tmp,'probabilities'),1,max)
}

# identify wrong classified trn.data~bag pair
wrong.classify<-which(prediction!=matrix(rep(trn.labels,B),ncol=B,byrow=FALSE),arr.ind=TRUE)
colnames(wrong.classify)<-c('data','bag')
wrong.classify.mat<-array(TRUE,c(nrow(trn.data),B))
wrong.classify.mat[wrong.classify]=FALSE

# self pruning
self.pruning.mat<-array(NA,c(nrow(trn.data),B))
for(i in 1:B)
	self.pruning.mat[,i]<-is.element(1:nrow(trn.data),idx.boot[,i])

self.pruning<-which(self.pruning.mat==TRUE,arr.ind=TRUE)
colnames(self.pruning)<-c('data','bag')

# make probability at wrong classify infinity
probability[wrong.classify]=-Inf
# make probability at self contained infinity
probability[self.pruning]=-Inf


# sort probability
r=21
probability.idx<-probability.dist<-array(NA,c(nrow(trn.data),r))
for (i in 1:nrow(trn.data)){
	probability.idx[i,]<-order(probability[i,],decreasing=TRUE)[1:r]
	probability.dist[i,]<-probability[i,probability.idx[i,]]
}



# knn pruning
k=31
knn.dist.mat<-my_knn_dist(trn.data,tst.data,k=k) #|tst|*|trn|(k)

# testing
tst.predict<-array(NA,c(nrow(tst.data),k*r))
for (i in 1:nrow(tst.data)){
	ith.trn.idx<-knn.dist.mat$knn.idx[i,]
	ith.bag.idx<-probability.idx[ith.trn.idx,]
	counter=1
	for (b in ith.bag.idx){	
		tmp<-predict(model[[b]],tst.data)
		tst.predict[,counter]<-as.numeric(levels(tmp)[tmp])
		counter=counter+1
	}
}

#aggregate
aggregated.tst.predict<-my_majority_vote(tst.predict[,1:k*r])
#error measure
my_classification_measure(aggregated.tst.predict$value,tst.labels)


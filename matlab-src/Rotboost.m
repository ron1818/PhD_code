function [Y_pred,Y_votes]=Rotboost(TrainX,TrainY,TestX,TestY,size1,size2)
% size1 number of rotations
% size2 number of trees for each rotation
   model=cell(size1*size2,1);
   Y_pred=zeros(length(TestY),1);
   Y_votes=zeros(length(TestY),size1*size2);
   alfa=zeros(1,size2*size1);
    ratio=0.75;
    K=3;
    testXRFnew=cell(size1,1);
 for loop=1:size1
     fprintf([num2str(loop),'th rotation\n'])
   [~,~,trainRFnew,testXRFnew{loop}]=RotationFal(TrainX, TrainY, TestX, K, ratio);
   [m,~]=size(trainRFnew);
    D=zeros(m,size2+1);
    D(:,1)=1/m;
    err=zeros(1,size2);
   
   i=1;

while (i<=size2)
    
    fprintf([num2str(i),'th decision in ', num2str(loop),'th rotatio\n'])
   
   index=randsample(1:m,m,true,D(:,i)); 
    model{(loop-1)*size2+i}=Obliquecartree_train(trainRFnew(index,:),TrainY(index),index,'nvartosample',round(sqrt(size(trainRFnew,2))),'oblique',4);
    %  Obliquecartree_train(Data,Labels,Index,varargin)
  Y_temp=Obliquecartree_predict(trainRFnew,trainRFnew,model{(loop-1)*size2+i});
  Y_votes(:,(loop-1)*size2+i)=Obliquecartree_predict(trainRFnew,testXRFnew{loop},model{(loop-1)*size2+i});
 
   err(i)=length(find(Y_temp'~=TrainY))/length(TrainY);
   
   if err(i)<0.5
    if err(i)==0
        err(i)=exp(-10);
    end
    alfa((loop-1)*size2+i)=0.5*log((1-err(i))/err(i));
    index_pos=find(Y_temp'==TrainY);
    index_neg=find(Y_temp'~=TrainY);
    D_temp=D(:,i);
    D_temp(index_pos)=D_temp(index_pos)*exp(-alfa(i));
    D_temp(index_neg)=D_temp(index_neg)*exp(alfa(i));
    D(:,i+1)=D_temp/sum(D_temp);
    i=i+1;
   else
       D(:,i)=1/m;
        
    end
end
 
unique_label=unique(TrainY);
nClass=length(unique_label);
for j=1:length(TestY)
    
    votes=zeros(size1,1);
    for jj=1:size1
        Y_temp2=Y_votes(j,(jj-1)*size2+1:(jj)*size2);
        alfa_temp=alfa((jj-1)*size2+1:(jj)*size2);
        counts=zeros(1,nClass);
        for j_j=1:nClass
           index= Y_temp2==unique_label(j_j) ;
           
           counts(j_j)=sum(alfa_temp(index));
        end
        [~,max_index]=max(counts);
       votes(max_index)=votes(max_index)+1;
    end
    [~,max_index1]=max(votes);
    Y_pred(j)=max_index1;
end
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
    end
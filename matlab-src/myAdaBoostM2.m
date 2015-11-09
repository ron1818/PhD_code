function [ z_pred, z_table, beta ] = myAdaBoostM2( x, y, z, maxIter )
% Experiments with a New Boosting Algorithm
% Yoav Freund Robert E. Schapire

% example 1
% load fisheriris
% m=size(meas,1);
% trnIdx=randsample(m,round(m/2));
% tstIdx=setdiff((1:m)',trnIdx);
% x=meas(trnIdx,:);
% y=species(trnIdx,:);
% z=meas(tstIdx,:);
% z_obs=species(tstIdx,:);
% maxIter=50;
% [ z_pred, z_table, beta ] = myAdaBoostM2( x, y, z, maxIter );

% example 1
% load ionosphere
% m=size(meas,1);
% trnIdx=randsample(m,round(m/2));
% tstIdx=setdiff((1:m)',trnIdx);
% x=meas(trnIdx,:);
% y=species(trnIdx,:);
% z=meas(tstIdx,:);
% z_obs=species(tstIdx,:);
% maxIter=50;
% [ z_pred, z_table, beta ] = myAdaBoostM2( x, y, z, maxIter );

y=nominal(y); % convert to nominal type
categories_y=categories(y); % categories of y
k=size(categories_y,1); % number of classes
[m,n]=size(x);

B=repmat(y,1,k)~=repmat(nominal(categories_y)',m,1); % 0 means yi (correct), 1 means y (incorrect)
y_map=not(B);
% initialize weight distr and error
D=repmat(1/m, m, 1); 
% weight vector, w_i,y
w_new=B.*repmat(D./(k-1),1,k);

T=maxIter;

e=zeros(1,T);
beta=ones(1,T);
cart=cell(1,T);
% iterate
t=1;
th_e=0;
while t<=T
    w=w_new;
    % W_i=\sum_y\ne y_i w_i,y
    W=sum(w, 2);
    % q(i,y)=w_y,y/W_i
    q=w./repmat(W,1,k);
    % D(i)=W/\sum W
    D=W./sum(W);
    % bootstrap with dist D
    [tmp,bootsam]=bootstrp(1, @mean, x, 'Weights',D); 
    tth_x=x(bootsam,:);
    tth_y=y(bootsam,:);
    % call weak learner, CART
    cart{t}=classregtree(tth_x,tth_y,'minparent', m/2);
    % getback to hypothesis
    y_fit=eval(cart{t},x);
    y_fit=nominal(y_fit);

    incorrect_idx=y_fit~=y;% incorrect prediction
    y_fit_map=repmat(y_fit,1,k)==repmat(nominal(categories_y)',m,1);
    incorrect_pred=repmat(incorrect_idx,1,k).*y_fit_map;
    incorrect_weight=sum(q.*incorrect_pred,2)+1;
    incorrect_weight=incorrect_idx.*incorrect_weight;
    e(t)=0.5*D'*incorrect_weight;
    if(e(t)==0) % perfect bypass this one
        continue;
    end
    % calculate error rate
    th_e=e(t);
    % weight updating param
    beta(t)=e(t)/(1-e(t));
    % update w
    beta_mat=y_fit_map-y_map;
    one_idx=beta_mat==1;
    zero_idx=beta_mat==0;
    negone_idx=beta_mat==-1;
    beta_mat(one_idx)=0;
    beta_mat(zero_idx)=1;
    beta_mat(negone_idx)=0.5;
    w_new=w.*beta(t).^beta_mat;
    t=t+1;
end
max_t=t-1;
% final model: cart, beta(t)
ww=log(1./beta(1:max_t));

% decision table
z_table=nominal(zeros(size(z,1),max_t)); % k by T
for t=1:max_t
    z_pred=eval(cart{t},z);
    z_pred=nominal(z_pred);
    z_table(:,t)=z_pred;
end
numeric_table=zeros(k,size(z_table,2));
z_pred=nominal(zeros(size(z,1),1));
for i=1:size(z_table,1) % each sample
    for t=1:max_t % each prdictor
        numeric_table(:,t)=z_table(i,t)==categories_y;
    end
    fin_table=repmat(ww,k,1).*numeric_table;
    z_fin=sum(fin_table,2);
    z_pred(i)=categories_y(z_fin==max(z_fin));
    z_pred(i)=nominal(z_pred(i));
end
end
        


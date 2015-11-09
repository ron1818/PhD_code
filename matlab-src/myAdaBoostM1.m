function [ z_pred, z_table, beta ] = myAdaBoostM1( x, y, z, maxIter )
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
% [ z_pred, z_table, beta ] = myAdaBoostM1( x, y, z, maxIter );

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
% [ z_pred, z_table, beta ] = myAdaBoostM1( x, y, z, maxIter );


y=nominal(y); % convert to nominal type
categories_y=categories(y); % categories of y
k=size(categories_y,1); % number of classes
[m,n]=size(x);
%initialize weight distr and error
D_new=repmat(1/m, m, 1);

T=maxIter; % max iter
eth=0.5; % error threshold

e=zeros(1,T);
beta=ones(1,T);
cart=cell(1,T);
% iterate
t=1;
th_e=0;
while t<=T && th_e<eth
    D=D_new;
    % bootstrap with dist D
    [tmp,bootsam]=bootstrp(1, @mean, x, 'Weights',D);
    tth_x=x(bootsam,:);
    tth_y=y(bootsam,:);
    % call weak learner, CART
    cart{t}=classregtree(tth_x,tth_y, 'minparent', m/2);
    % getback to hypothesis
    y_fit=eval(cart{t},x);
    y_fit=nominal(y_fit);
    % calculate error rate
    e(t)=D'*(y_fit~=y);
    if(e(t)==0) % perfect bypass this one
        continue;
    end
    % question!!!
    % if all correct (perfect predict), then e=0 and beta=0
    th_e=e(t);
    % weight updating param
    beta(t)=e(t)/(1-e(t));
    % update distr
    correct=beta(t)*(y_fit==y);
    incorrect=(y_fit~=y);
    D_new=D.*(correct+incorrect);
    % normalize
    D_new=D_new./sum(D_new);
    t=t+1;
end
max_t=t-1;
% final model: cart, beta(t)
ww=log(1./beta(1:max_t));


% decision table
z_table=nominal(zeros(size(z,1),max_t)); % k by T
for t=1:max_t
    pred=eval(cart{t},z);
    pred=nominal(pred);
    z_table(:,t)=pred;
end
numeric_table=zeros(k,size(z_table,2));
z_pred=nominal(zeros(size(z,1),1));
for i=1:size(z_table,1) % each sample
    for t=1:size(z_table,2) % each prdictor
        numeric_table(:,t)=z_table(i,t)==categories_y;
    end
    fin_table=repmat(ww,k,1).*numeric_table;
    z_fin=sum(fin_table,2);
    tmp=categories_y(z_fin==max(z_fin));
    z_pred(i)=nominal(tmp);
end
end

function [ z_pred, beta, f ] = myAdaBoostR2_ANN( x, y, z, maxIter, pcoef, hiddenLayerSize, K )
% Experiments with a New Boosting Algorithm
% Yoav Freund Robert E. Schapire

% use ANN as predict algo

% % example 1
% load moore
% y = moore(:,6);              % Response
% x = moore(:,1:5);            % Original predictors
% z = x+10*randn(size(x));       % Correlated predictors
% pcoef=2;
% maxIter=50;

% % example 2
% load acetylene
% % scale data
% [s_x1] = scale_data( x1,1,0,[],[] );
% [s_x2] = scale_data( x2,1,0,[],[] );
% [s_x3] = scale_data( x3,1,0,[],[] );
% [s_y] = scale_data( y,1,0,[],[] );
% trnIdx=randsample(16,round(16/3));
% tstIdx=setdiff((1:16)',trnIdx);
% Y = s_y(trnIdx,:);              % Response
% X = [s_x1(trnIdx,:), s_x2(trnIdx,:), s_x3(trnIdx,:)]; % Original predictors
% z = [s_x1(tstIdx,:), s_x2(tstIdx,:), s_x3(tstIdx,:)]; 
% z_obs=s_y(tstIdx,:);
% pcoef=2;
% maxIter=50;
% [ z_pred, beta ] = myAdaBoostR2( X, Y, z, maxIter, pcoef )


[m,n]=size(x);
%initialize weight distr and error
D_new=repmat(1/m, m, 1);

T=maxIter; % max iter
L_bar_new=0; % average loss fn

beta=ones(1,T);
f=cell(1,T);
% iterate
t=1;
while t<=T && L_bar_new<0.5
    D=D_new;
    L_bar=L_bar_new;
    % bootstrap with dist D
    [tmp,bootsam]=bootstrp(1, @mean, x, 'Weights',D);
    tth_x=x(bootsam,:);
    tth_y=y(bootsam,:);
    % call myANN
     [ y_fit, f{t}, best_CV_MSE ] = myANN( tth_x, tth_y, x, y, hiddenLayerSize, 'trainlm', 'tansig', K);
%     % call weak learner, CART
%     f{t}=trn_fun(tth_x,tth_y);
%     % getback to hypothesis
%     y_fit=tst_fun(f{t},x);
    l=abs(y_fit-y); %loss for each training example
    denorm=max(l);
    % Calculate the loss function for each training example using three different functional forms
    switch pcoef
        case 1, % linear
            L=l/denorm;
        case 2, %square law
            L=(l/denorm).^2;
        case 3, % exponential
            L=1-exp(-l/denorm);
        otherwise, % default linear
            L=l/denorm;
    end
    
    L_bar_new=L'*D; % update average loss
    beta(t)=L_bar_new/(1-L_bar_new);
    D_new=D.*beta(t).^(1-L);    
    % normalize
    D_new=D_new./sum(D_new);
    
    t=t+1;
end
max_t=t-1;
% final model: cart, beta(t)
ww=log(1./beta(1:max_t));

% test
pred=zeros(size(z,1), max_t);
for t=1:max_t
    pred(:,t)=f{t}(z')';
end

% weighted median using ww as final hypo
weighted_pred=repmat(ww, size(z,1),1).*pred;

z_pred=median(weighted_pred,2)/median(ww);
end




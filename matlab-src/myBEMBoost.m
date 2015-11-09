function [ z_pred ] = myBEMBoost( trn_fun, tst_fun, x, y, z, maxIter, BEM )
% Feely 2000

[m,n]=size(x);
%initialize weight distr and error
D_new=repmat(1/m, m, 1);

T=maxIter; % max iter
% phi=0.5; % error threshold

errorcnt=zeros(1,T);
beta=ones(1,T);
f=cell(1,T);
% iterate
t=1;
while t<=T
    D=D_new;
    % bootstrap with dist D
    [tmp,bootsam]=bootstrp(1, @mean, x, 'Weights',D);
    tth_x=x(bootsam,:);
    tth_y=y(bootsam,:);
    % call weak learner, CART
    f{t}=trn_fun(tth_x,tth_y);
    % getback to hypothesis
    y_fit=tst_fun(f{t},x);
    AE=abs(y_fit-y); % absolute error
    % calculate error count AE>BEM
    is_ae_lt_BEM=AE>BEM;
    errorcnt(t)=sum(is_ae_lt_BEM);
    % calculate upfactor and downfactor
    upfactor(t)=m/errorcnt(t);
    downfactor(t)=1/upfactor(t);
%     
%     if(e(t)==0) % perfect bypass this one
%         continue;
%     end
    
    % update distr
    is_ae_le_BEM=not(is_ae_lt_BEM);
    
    D_new=D.*(is_ae_lt_BEM.*upfactor(t)+is_ae_le_BEM.*downfactor(t));
    % normalize
    D_new=D_new./sum(D_new);
    t=t+1;
end
max_t=t-1;
% final model: cart, beta(t)
% ww=log(1./beta(1:max_t));

% test
pred=zeros(size(z,1), max_t);
for t=1:max_t
    pred(:,t)=tst_fun(f{t},z);
end


z_pred=median(pred,2);
end




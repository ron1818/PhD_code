function [ z_pred, beta ] = myAdaBoostRT( trn_fun, tst_fun, x, y, z, maxIter, phi, pcoef )
% Experiments with a New Boosting Algorithm
% Yoav Freund Robert E. Schapire
% pcoef=2;

% example 1
% load moore
% y = moore(:,6);              % Response
% x = moore(:,1:5);            % Original predictors
% z = x+10*randn(size(x));       % Correlated predictors
% phi=0.5;
% maxIter=50;
% pcoef=2;
% [ z_pred, beta ] = myAdaBoostRT( x, y, z, maxIter, phi, pcoef );

[m,n]=size(x);
%initialize weight distr and error
D_new=repmat(1/m, m, 1);

T=maxIter; % max iter
% phi=0.5; % error threshold

e=zeros(1,T);
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
    ARE=abs((y_fit-y)./y);
    is_are_lt_phi=ARE>phi;
    
    % calculate error rate
    e(t)=D'*is_are_lt_phi;
    if(e(t)==0) % perfect bypass this one
        continue;
    end
    % question!!!
    % if all correct (perfect predict), then e=0 and beta=0
    %     th_e=e(t);
    % weight updating param
    beta(t)=e(t).^pcoef;
    % update distr
    lt_phi=beta(t)*(ARE<=phi);
    ge_phi=is_are_lt_phi;
    D_new=D.*(lt_phi+ge_phi);
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
    pred(:,t)=tst_fun(f{t},z);
end

norminator=ww*pred';
denorminator=sum(ww);

z_pred=norminator./denorminator;
z_pred=z_pred';
end




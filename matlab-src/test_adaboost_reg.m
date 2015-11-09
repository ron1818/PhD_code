% addpath
addpath('../libsvm-3.11/');
%%
load moore
org_Y = moore(:,6);              % Response
org_X = moore(:,1:5);            % Original predictors
X=scale_data( org_X,1,0,[],[] );
Y=scale_data( org_Y,1,0.01,[],[] );
z = X+0.1*randn(size(X));       % Correlated predictors

ens=fitensemble(X,Y,'LSBoost',100,'Tree');
ens_z=ens.predict(z);

pcoef=1;
maxIter=100;
[ z_pred_r2, beta_r2 ] = myAdaBoostR2( @classregtree, @eval, X, Y, z, maxIter, pcoef );
phi=0.5;
[ z_pred_rt, beta_rt ] = myAdaBoostRT( @classregtree, @eval, X, Y, z, maxIter, phi, pcoef );
BEM=0.1;
[ z_pred_bem ] = myBEMBoost( @classregtree, @eval, X, Y, z, maxIter, BEM );

plot([z_pred_r2, z_pred_rt, z_pred_bem, ens_z, Y]);
legend('R2', 'RT', 'BEM', 'toolbox', 'Obs');


%%
load acetylene
% scale data
[s_x1] = scale_data( x1,1,0,[],[] );
[s_x2] = scale_data( x2,1,0,[],[] );
[s_x3] = scale_data( x3,1,0,[],[] );
[s_y] = scale_data( y,1,0.01,[],[] );
trnIdx=randsample(16,2*round(16/3));
tstIdx=setdiff((1:16)',trnIdx);
Y = s_y(trnIdx,:);              % Response
X = [s_x1(trnIdx,:), s_x2(trnIdx,:), s_x3(trnIdx,:)]; % Original predictors
z = [s_x1(tstIdx,:), s_x2(tstIdx,:), s_x3(tstIdx,:)]; 
z_obs=s_y(tstIdx,:);

ens=fitensemble(X,Y,'LSBoost',100,'Tree');
ens_z=ens.predict(z);

pcoef=1;
maxIter=100;
[ z_pred_r2, beta_r2 ] = myAdaBoostR2( @classregtree, X, Y, z, maxIter, pcoef );
phi=0.5;
[ z_pred_rt, beta_rt ] = myAdaBoostRT( @classregtree, X, Y, z, maxIter, phi, pcoef );
BEM=0.1;
[ z_pred_bem ] = myBEMBoost( @classregtree, X, Y, z, maxIter, BEM );

plot([z_pred_r2, z_pred_rt, z_pred_bem, ens_z, Y]);
legend('R2', 'RT', 'BEM', 'toolbox', 'Obs');

%%
load carbig
org_x=[Horsepower Weight Displacement Acceleration Cylinders];
org_y=MPG;
% remove nan
nan_idx=isnan(org_y);
org_x(nan_idx,:)=[];
org_y(nan_idx,:)=[];
% scale data
x=scale_data( org_x,1,0.01,[],[] );
y=scale_data( org_y,1,0.01,[],[] );
trnIdx=randsample((406-length(nan_idx)),200);
tstIdx=setdiff((1:(406-length(nan_idx)))',trnIdx);
X=x(trnIdx,:);
Y=y(trnIdx,:);
z =x(tstIdx,:);
z_obs=y(tstIdx,:);

ens=fitensemble(X,Y,'LSBoost',100,'Tree');
ens_z=ens.predict(z);

pcoef=1;
maxIter=100;
[ z_pred_r2, beta_r2 ] = myAdaBoostR2( @classregtree, X, Y, z, maxIter, pcoef );
phi=0.5;
[ z_pred_rt, beta_rt ] = myAdaBoostRT( @classregtree, X, Y, z, maxIter, phi, pcoef );
BEM=0.1;
[ z_pred_bem ] = myBEMBoost( @classregtree, X, Y, z, maxIter, BEM );

plot([z_pred_r2, z_pred_rt, z_pred_bem, ens_z, Y]);
legend('R2', 'RT', 'BEM', 'toolbox', 'Obs');




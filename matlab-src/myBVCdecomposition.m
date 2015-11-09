%% BV decomp
% bias: amount we are off the real function
% variance: variance of the predictions
% E(square error)=var{noise}+bias^2+var{pred}
% E[(t_i-y_i)]=E[\epsilon^2]+E[(f_i-E[y_i])^2]+E[(E[y_i]-y_i)^2]

y_pred=rand(100,10); % prediction
y=rand(100,1); % observation

mse=mean((y-y_pred).^2); % mse of each predictor

bias=mean((y-mean(y_pred)).^2);
variance=mean((y_pred-mean(y_pred)).^2);
error_variance=mse-bias^2-variance;

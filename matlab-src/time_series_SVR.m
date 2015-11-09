% this script requires a 1 dimensional time series and SVR functions
% svmtrain and svmpredict.
% the time series used here is a random vector

x_trn=rand(100,1); % training time series
x_tst=rand(20,1); % testing time series

HORIZON=1; % one step ahead
F_COUNT=4; % vector length
ispoint=1; % trn_labels is single value
[ trn_data, trn_labels ]=ts2mat(x_trn, HORIZON, F_COUNT, ispoint);
[ tst_data, tst_labels ]=ts2mat(x_tst, HORIZON, F_COUNT, ispoint);

% SVR train
options=sprintf('-s 3 -t 2');
model=svmtrain(trn_labels, trn_data, options);
% SVR predict
[tst_pred, accuracy, decision]=svmpredict(tst_labels, tst_data, model);
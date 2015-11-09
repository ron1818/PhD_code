function [ aggregated_predicted_labels, predicted_labels ]= bagging_ANN( X_trn, Y_trn, X_tst, Y_tst, hiddenLayerSize, bag_numbers )
% bootstraping
if license('checkout', 'Statistics_Toolbox')
    [tmp, bags] = bootstrp(bag_numbers, [], Y_trn, X_trn);
elseif exist('randi', 'builtin')
    bags=randi(length(Y_trn), length(Y_trn), bag_numbers);
else
    bags=randint(length(Y_trn), bag_numbers, [1 length(Y_trn)]);
end
clear tmp;

% bagging training
% preallocation
predicted_labels=zeros(length(Y_tst), bag_numbers);
OOB=zeros(1, bag_numbers);
OOB_idx=zeros(length(Y_trn), bag_numbers);
for i=1:bag_numbers % each bag
    ith_bag=bags(:,i);
    ith_data=X_trn(ith_bag,:);
    ith_labels=Y_trn(ith_bag);
    ith_outbag=setdiff((1:length(Y_trn)),ith_bag);
    ith_outbag_data=X_trn(ith_outbag,:);
    OOB_idx(ith_outbag,i)=1;
    ith_outbag_labels=Y_trn(ith_outbag);
    ith_predicted_OOB=myANN( ith_data, ith_labels, ith_outbag_data, ith_outbag_labels, hiddenLayerSize);
    OOB(i)=sqrt(mean((ith_predicted_OOB - ith_outbag_labels).^2)); % RMSE
    predicted_labels(:,i)=myANN( ith_data, ith_labels, X_tst, Y_tst, hiddenLayerSize);
end
% aggregate
aggregated_predicted_labels=mean(predicted_labels,2);
end




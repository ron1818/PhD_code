function[ trn_data, trn_labels, tst_data, tst_labels ]=ts2mat(vector, split_point, HORIZON, F_COUNT, ispoint)
% UPDATED 20140327
% if split_point==0, no split, all go to trn_data/labels, no
% tst_data/labels
%
% UPDATED 20140110
% change F_COUNT to a vector so that discontinuous vector can be used in
% trn_data and testing data
%
% UPDATED 20131226
% using sliding window for ts2mat transformation, no need separation of
% training and testing TS.
% input: split point: if <1 (fraction), then it means percentage
% if >1 (integer), then it means spliting index
% outputs: trn_data, trn_labels, tst_data, tst_labels
%
% convert ts to matrix with h lags
% inputs:
% vector: n*1 ts
% HORIZON: forecasting horizon
% F_COUNT: matrix col number
% ispoint: point forecasting labels=m*1 vector, otherwise, labels=m*H vector
% outputs:
% matrix: ts in matrix format: m*F
% labels: m*1 or m*H
% example:
% vector=(1:100)';
% HORIZON=2;
% F_COUNT=4;
% ispoint=1;
% [ matrix, labels ]=ts2mat(vector, HORIZON, F_COUNT, ispoint);

max_HORIZON=max(HORIZON);

length_F_COUNT=length(F_COUNT);
max_F_COUNT=max(F_COUNT);
matrix_length=length(vector)-max_F_COUNT-max_HORIZON+1;
matrix=zeros(matrix_length, length_F_COUNT);

labels=zeros(matrix_length, max_HORIZON);

for i =1:max_HORIZON
    labels(:,i)=vector((i+max_F_COUNT):(i+max_F_COUNT+matrix_length-1));
end

for i =1:length_F_COUNT
    matrix(:,i)=vector(F_COUNT(i):(F_COUNT(i)+matrix_length-1));
end

% split to trn and tst
if split_point==0 % no split
    trn_length=size(matrix,1);
elseif split_point<1 % fraction
    trn_length=round(split_point*length(vector))-max_F_COUNT;
else % integer
    trn_length=split_point;
end
trn_data=matrix(1:trn_length,:);
trn_labels=labels(1:trn_length,:);

if split_point==0;
    tst_data=[];
    tst_labels=[];
else
    tst_data=matrix((trn_length+1):end,:);
    tst_labels=labels((trn_length+1):end,:);
end

if ispoint
    trn_labels=trn_labels(:,HORIZON);
    tst_labels=tst_labels(:,HORIZON);
end
end


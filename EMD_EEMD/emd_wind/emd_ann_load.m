%%
addpath('../dataset/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
% result_date=date;
% HORIZON_array=[1, 2, 12, 24, 48]; % 0.5, 1, 6, 12, 24 hour
HORIZON_array=[48]; % 24 hour
F_COUNT=[1:3, 23:25, 46:48];
RATIO=[0.7, 0.3];

% ANN_predict=cell(1, length(week_idx));
% KNN_best_k=zeros(length(week_idx), max(HORIZON_array));
% KNN_best_d=zeros(length(week_idx), max(HORIZON_array));
% KNN_mat=cell(1, max(HORIZON_array));
%
% EMDKNN_predict=cell(1, length(week_idx));
% EMDKNN_best_k=cell(length(week_idx),max(HORIZON_array));
% EMDKNN_best_d=cell(length(week_idx),max(HORIZON_array));
% EMDKNN_mat=cell(6, max(HORIZON_array));
%
% EAKNN_predict=cell(1, length(week_idx));
% EAKNN_best_k=zeros(length(week_idx), max(HORIZON_array));
% EAKNN_best_d=zeros(length(week_idx), max(HORIZON_array));
% EAKNN_mat=cell(1, max(HORIZON_array));

% result=cell(max(week_idx),1);

import_data=importdata('countryenergy_NSW_AU_2013.csv'); % read data
load_data=import_data.data; % discard timestamp
load_data=na_correction(load_data); % correct na entry

% skip outlier detection
data=load_data;

%     %% detect outliers
%     IQR_factor=1.5;
%     window=30;
%     threshold=10;
%     outlier_idx1=residual_IQR_outlier( load_data, IQR_factor );
%     outlier_idx2=window_mad_outlier( load_data, window, threshold );
%     % remove outliers
%     outlier_idx=intersect(outlier_idx1, outlier_idx2);
%     data = outlier_correction( load_data, outlier_idx );

%% EMD
% scale
[scaled_data, max_data, min_data]=scale_data( data,1,0,[],[] );

% EMD IMF series
options.MAXMODES=5;
IMF=emd(scaled_data, options);
IMF=IMF';

%% split training and testing
ispoint=0;
% split original ts
[ trn_data, trn_labels, tst_data, tst_labels ]=ts2mat(scaled_data, RATIO(1), max(HORIZON_array), F_COUNT, ispoint);

% %split IMF
% [ trn_data, trn_labels, tst_data, tst_labels ]=ts2mat(data, RATIO(1), max(HORIZON_array), F_COUNT, ispoint);



%     ha = tight_subplot(6, 1, [0.01 0.01], [.1 .01],[.1 .01]);
%     for k=1:6
%         plot(IMF_trn(:,k));
%         ylim(ha(k), [min(IMF_trn(:,k)), max(IMF_trn(:,k))]);
%     end


% knn parameter
%     k_range=2*round(sqrt(length(trn_speed)-12)/2)+(-11:2:11);
%     k_range=abs(k_range);
%     d_range=4:2:9;
%     w_flag=1;


%% Persistent
per_residue=tst_labels-repmat(tst_data(:,end), 1, max(HORIZON_array));
[persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(tst_data(:,end), 1, max(HORIZON_array)), per_residue);
persistent_RMSE=persistent_RMSE(HORIZON_array);
persistent_sMAPE=persistent_sMAPE(HORIZON_array);
persistent_MASE=persistent_MASE(HORIZON_array);

%% ANN
ANN_RMSE=zeros(1,size(HORIZON_array,2));
ANN_sMAPE=zeros(1,size(HORIZON_array,2));
ANN_MASE=zeros(1,size(HORIZON_array,2));
counter=1;
for HORIZON=HORIZON_array
    hiddenlayer=2*size(trn_data,2);
    [ANN_predict]=myANN(trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON),hiddenlayer); % hidden layer =10
    [ANN_RMSE(counter), ANN_sMAPE(counter), ANN_MASE(counter)]=myErrorMeasure(tst_labels(:,HORIZON), ANN_predict, per_residue(:,HORIZON));
    counter=counter+1;
end



%% EMD-ANN
% EMDANN_predict=zeros(length(tst_labels), max(HORIZON_array));
for HORIZON=HORIZON_array
%     ith_predict=zeros(length(tst_labels), size(IMF,1));
    for i=1:size(IMF,2)
        i
        [ith_trn_data, ith_trn_labels, ith_tst_data, ith_tst_labels] = ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
        hiddenlayer=2*size(ith_trn_data,2);
        [tmp_predict]=myANN(ith_trn_data, ith_trn_labels(:,HORIZON), ith_tst_data, ith_tst_labels(:,HORIZON),hiddenlayer); % hidden layer =10
        ith_predict(:,i)=tmp_predict;
%         ith_predict(:,i)=tmp_predict(length(tmp_predict)-length(ith_tst_labels)+1:end);
    end
    EMDANN_predict=sum(ith_predict,2);
end

[EMD_ANN_RMSE, EMD_ANN_sMAPE, EMD_ANN_MASE]=myErrorMeasure(tst_labels, EMDANN_predict, per_residue);
EMD_ANN_RMSE=EMD_ANN_RMSE(HORIZON_array);
EMD_ANN_sMAPE=EMD_ANN_sMAPE(HORIZON_array);
EMD_ANN_MASE=EMD_ANN_MASE(HORIZON_array);

%% EMD-A-ANN

% for HORIZON=HORIZON_array
    ith_predict=zeros(length(tst_labels), size(IMF,1));
    for i=1:size(IMF,2)
        step=length(F_COUNT);
        [imf_trn_data(:,(step*(i-1)+(1:step))), imf_trn_labels, imf_tst_data(:,(step*(i-1)+(1:step))), imf_tst_labels] = ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
    end
    counter=1;
    for HORIZON=HORIZON_array
        hiddenlayer=2*size(imf_trn_data,2);
        [EMDAANN_predict(:,counter)]=myANN(imf_trn_data, trn_labels(:,HORIZON), imf_tst_data, tst_labels(:,HORIZON),50); % hidden layer =10
        counter=counter+1;
    end

[EMD_AANN_RMSE, EMD_AANN_sMAPE, EMD_AANN_MASE]=myErrorMeasure(tst_labels(:,HORIZON_array), EMDAANN_predict, per_residue(:,HORIZON_array));

%%
result=[persistent_RMSE, persistent_sMAPE, persistent_MASE
    ANN_RMSE, ANN_sMAPE, ANN_MASE
    EMD_ANN_RMSE, EMD_ANN_sMAPE, EMD_ANN_MASE
    EMD_AANN_RMSE, EMD_AANN_sMAPE, EMD_AANN_MASE];

save('result_emdann_20140115.mat', 'result', 'ANN_predict', 'EMDANN_predict', 'EMDAANN_predict');
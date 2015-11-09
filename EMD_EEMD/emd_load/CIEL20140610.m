%%
addpath('../dataset/load/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
HORIZON_array=[12 24 48];
ispoint=0;
numLags=96;
% import data
import_data=importdata('countryenergynsw_2011.csv'); % read data
time=datestr(import_data.textdata(2:end,1));

% correct na entry
load_data=na_correction(import_data.data);

%   skip
%   %% outlier detection
%   data=load_data;
%   skip
%   %% detect outliers
%   IQR_factor=1.5;
%   window=30;
%   threshold=10;
%   outlier_idx1=residual_IQR_outlier( load_data, IQR_factor );
%   outlier_idx2=window_mad_outlier( load_data, window, threshold );
%   % remove outliers
%   outlier_idx=intersect(outlier_idx1, outlier_idx2);
%   data = outlier_correction( load_data, outlier_idx );

% scale
[total_scaled_data, max_data, min_data]=scale_data( load_data,1,0,[],[] );
ts=timeseries(total_scaled_data, time);

% prepare for 12 typical months
starttime={'2011-01-01 00:00', '2011-02-01 00:00', '2011-03-01 00:00',...
           '2011-04-01 00:00', '2011-05-01 00:00', '2011-06-01 00:00', '2011-07-01 00:00',...
           '2011-08-01 00:00' '2011-09-01 00:00' '2011-10-01 00:00' '2011-11-01 00:00', '2011-12-01 00:00'};
endtime={'2011-01-31 23:30', '2011-02-28 23:30', '2011-03-31 23:30',...
           '2011-04-30 23:30', '2011-05-31 23:30', '2011-06-30 23:30', '2011-07-31 23:30',...
           '2011-08-31 23:30' '2011-09-30 23:30' '2011-10-31 23:30' '2011-11-30 23:50', '2011-12-31 23:30'};

% preallocate
result=cell(1,12);
pred=cell(1,12);
delays=cell(1,12);
SVR_param=cell(1,12);
for idx=1:12 % 12 typical months, 6 only for conference
    % form a TS
    scaled_ts=getsampleusingtime(ts, starttime{idx}, endtime{idx});
    scaled_data=scaled_ts.data;
    
    % use pacf to determine the feedback Delays.
    par_trn=parcorr(scaled_data, numLags);
    % select the top 10% quantile
    feedbackDelays=find(abs(par_trn)>quantile(abs(par_trn), 0.9));
    % if use pacf threshold:
    % feedbackDelays=find(abs(par_trn)>1.96/sqrt(length(trnIdx)));
    feedbackDelays=feedbackDelays-1;
    feedbackDelays=feedbackDelays(2:end);
    hiddenLayerSize=2*length(feedbackDelays);
    
    RATIO=[24, 7]/31;
    
    [trn_data,trn_labels,tst_data,tst_labels]=ts2mat(scaled_data, (RATIO(1)), max(HORIZON_array), feedbackDelays,ispoint);
    TST_LB{idx}=tst_labels;
    trnIdx=1:(size(trn_data,1)+max(feedbackDelays));
    tstIdx=trnIdx(end)+(1:size(tst_data,1));
    
    % EMD IMF series
    IMF=emd(scaled_data,'MAXMODES',5);
    IMF=IMF';
    
    % pacf of imfs
    for i=1:size(IMF, 2)
        par=parcorr(IMF(1:round(size(IMF,1)*RATIO(1)),i),numLags);
        tmp_lag=find(abs(par)>1.96/sqrt(round(length(IMF(:,i))*RATIO(1))))-1;
        imf_lag{i}=setdiff(tmp_lag,0); % remove 0
    end

    
    for i=1:size(IMF, 2)
        [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), feedbackDelays,ispoint);
    end
  
    
    % plot IMF
%     myIMFplot(scaled_data(trnIdx), IMF(trnIdx,:));
    
    % statistical Signi test
%     myIMFstatistical_significance_test(IMF(trnIdx,:)', [0.05, 0.95]);
   
    %% Persistent
    % preallocation
    persistent_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    
    PERIOD=48;
    for HORIZON=HORIZON_array
        persistent_pred(:,HORIZON)=scaled_data(tstIdx-PERIOD+HORIZON-1);
    end
    per_residue=tst_labels-repmat(persistent_pred(:,end), 1, max(HORIZON_array));
    [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(persistent_pred(:,end), 1, max(HORIZON_array)), per_residue);
    
    
    %% BP
    % preallocation
    bpnn_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    
    for HORIZON=HORIZON_array
        [ bpnn_pred(:,HORIZON) ] = myANN( trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON), hiddenLayerSize);
    end
    [ BPNN_RMSE, BPNN_sMAPE, BPNN_MASE ]=myErrorMeasure(tst_labels, bpnn_pred, per_residue);
    
    %% AdaboostR2-CART
     % preallocation
    ab_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    maxIter=100;
    pcoef=2;
    for HORIZON=HORIZON_array
    [ ab_pred(:,HORIZON) ] = myAdaBoostR2( @classregtree, @eval, trn_data, trn_labels(:,HORIZON), tst_data, maxIter, pcoef );
    end
    [ AB_RMSE, AB_sMAPE, AB_MASE ]=myErrorMeasure(tst_labels, ab_pred, per_residue);
    
    %% AdaboostR2-BP
     % preallocation
    abpnn_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    maxIter=100;
    pcoef=2;
    K=2;
    for HORIZON=HORIZON_array
    [ abpnn_pred(:,HORIZON) ] = myAdaBoostR2_ANN( trn_data, trn_labels(:,HORIZON), tst_data, maxIter, pcoef, hiddenLayerSize, K );
    end
    [ ABPNN_RMSE, ABPNN_sMAPE, ABPNN_MASE ]=myErrorMeasure(tst_labels, abpnn_pred, per_residue);
    
    %% EMD-BP
    % preallocation
    ANN_IMF_pred=cell(1,size(IMF,2));
    for i=1:size(IMF, 2)
        % preallocation
        tmp_emd_ann_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        
        for HORIZON=HORIZON_array
            [ tmp_emd_ann_pred(:,HORIZON) ] = ...
                myANN( IMF_trn_data{i}, IMF_trn_labels{i}(:,HORIZON), IMF_tst_data{i}, IMF_tst_labels{i}(:,HORIZON), hiddenLayerSize );
        end
        ANN_IMF_pred{i}=tmp_emd_ann_pred;
    end
    mat_ANNIMF_pred=cell2mat(ANN_IMF_pred);
    
    % preallocation
    emdann_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    
    for i=1:max(HORIZON_array)
        [emdann_pred(:,i)]=sum(mat_ANNIMF_pred(:,i:max(HORIZON_array):end),2);
    end
    [ EMD_BPNN_RMSE, EMD_BPNN_sMAPE, EMD_BPNN_MASE ]=myErrorMeasure(tst_labels, emdann_pred, per_residue);
    
    %% EMD-AdaboostR2-ANN
    % preallocation
    ABPNN_IMF_pred=cell(1,size(IMF,2));
    maxIter=100;
    pcoef=2;
    K=2;
    for i=1:size(IMF, 2)
        % preallocation
        tmp_emd_abpnn_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        
        for HORIZON=HORIZON_array
            [ tmp_emd_abpnn_pred(:,HORIZON)]= myAdaBoostR2_ANN( IMF_trn_data{i}, IMF_trn_labels{i}(:,HORIZON), IMF_tst_data{i}, maxIter, pcoef, hiddenLayerSize, K );
        end
        ABPNN_IMF_pred{i}=tmp_emd_abpnn_pred;
    end
    mat_ABPNNIMF_pred=cell2mat(ABPNN_IMF_pred);
    
    % preallocation
    emdabpnn_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    
    for i=1:max(HORIZON_array)
        [emdabpnn_pred(:,i)]=sum(mat_ABPNNIMF_pred(:,i:max(HORIZON_array):end),2);
    end
    [ EMD_ABPNN_RMSE, EMD_ABPNN_sMAPE, EMD_ABPNN_MASE ]=myErrorMeasure(tst_labels, emdabpnn_pred, per_residue);
    
    %% EMD-AdaboostR2
    % preallocation
    AB_IMF_pred=cell(1,size(IMF,2));
    maxIter=100;
    pcoef=2;
    for i=1:size(IMF, 2)
        % preallocation
        tmp_emd_ab_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        
        for HORIZON=HORIZON_array
            [ tmp_emd_ab_pred(:,HORIZON)]= myAdaBoostR2( @classregtree, @eval, IMF_trn_data{i}, IMF_trn_labels{i}(:,HORIZON), IMF_tst_data{i}, maxIter, pcoef );
        end
        AB_IMF_pred{i}=tmp_emd_ab_pred;
    end
    mat_ABIMF_pred=cell2mat(AB_IMF_pred);
    
    % preallocation
    emdab_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    
    for i=1:max(HORIZON_array)
        [emdab_pred(:,i)]=sum(mat_ABIMF_pred(:,i:max(HORIZON_array):end),2);
    end
    [ EMD_AB_RMSE, EMD_AB_sMAPE, EMD_AB_MASE ]=myErrorMeasure(tst_labels, emdab_pred, per_residue);
    
   
    %%
    result{idx}=[persistent_RMSE, persistent_sMAPE, persistent_MASE
        AB_RMSE, AB_sMAPE, AB_MASE
        BPNN_RMSE, BPNN_sMAPE, BPNN_MASE
        ABPNN_RMSE, ABPNN_sMAPE, ABPNN_MASE 
        EMD_BPNN_RMSE, EMD_BPNN_sMAPE, EMD_BPNN_MASE
        EMD_AB_RMSE, EMD_AB_sMAPE, EMD_AB_MASE
        EMD_ABPNN_RMSE, EMD_ABPNN_sMAPE, EMD_ABPNN_MASE
        ];
    delays{idx}=feedbackDelays;
    pred{idx}=[persistent_pred, ab_pred, bpnn_pred, abpnn_pred, emdann_pred, emdab_pred, emdabpnn_pred];
end
save('result_CIEL_20140610.mat', 'result', 'TST_LB', 'delays', 'pred');

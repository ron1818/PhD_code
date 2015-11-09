%%
addpath('../dataset/load/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
HORIZON_array=[6 12 18 24 36 48];
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

% prepare for four typical months
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
for idx=1:12 % four typical months
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
    
    RATIO=[17, 7, 7]/31;
    
    [trn_data,trn_labels,tst_data,tst_labels]=ts2mat(scaled_data, (RATIO(1)+RATIO(2)), max(HORIZON_array), feedbackDelays,ispoint);
    
    trnIdx=1:(size(trn_data,1)+max(feedbackDelays));
    tstIdx=trnIdx(end)+(1:size(tst_data,1));
%     
%     % EMD IMF series
%     IMF=emd(scaled_data);
%     IMF=IMF';
%     
%     % CEEMD IMF series
%     Nstd=0.2;
%     NR=100;
%     MaxIter=500;
%     [ CEEMD_IMF ] = myCEEMD( scaled_data,Nstd,NR,MaxIter );
%     CEEMD_IMF=CEEMD_IMF';
%     
%     % CEEMDAN IMF series
%     NR=20;
%     [ CEEMDAN_IMF ] = ceemdan( scaled_data,Nstd,NR,MaxIter );
%     CEEMDAN_IMF=CEEMDAN_IMF';
%     
%     % pacf of imfs
%     for i=1:size(IMF, 2)
%         par=parcorr(IMF(1:round(size(IMF,1)*RATIO(1)),i),numLags);
%         tmp_lag=find(abs(par)>1.96/sqrt(round(length(IMF(:,i))*(RATIO(1)+RATIO(2)))))-1;
%         imf_lag{i}=setdiff(tmp_lag,0); % remove 0
%     end
%     
%     % pacf of imfs
%     for i=1:size(CEEMD_IMF, 2)
%         par=parcorr(CEEMD_IMF(1:round(size(CEEMD_IMF,1)*RATIO(1)),i),numLags);
%         tmp_lag=find(abs(par)>1.96/sqrt(round(length(CEEMD_IMF(:,i))*(RATIO(1)+RATIO(2)))))-1;
%         ceemd_imf_lag{i}=setdiff(tmp_lag,0); % remove 0
%     end
%     
%     % pacf of imfs
%     for i=1:size(CEEMDAN_IMF, 2)
%         par=parcorr(CEEMDAN_IMF(1:round(size(CEEMDAN_IMF,1)*RATIO(1)),i),numLags);
%         tmp_lag=find(abs(par)>1.96/sqrt(round(length(CEEMDAN_IMF(:,i))*(RATIO(1)+RATIO(2)))))-1;
%         ceemdan_imf_lag{i}=setdiff(tmp_lag,0); % remove 0
%     end
%     
%     for i=1:size(IMF, 2)
%         [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), (RATIO(1)+RATIO(2)), max(HORIZON_array), feedbackDelays,ispoint);
%     end
%     
%     for i=1:size(CEEMD_IMF, 2)
%         [ CEEMD_IMF_trn_data{i}, CEEMD_IMF_trn_labels{i}, CEEMD_IMF_tst_data{i}, CEEMD_IMF_tst_labels{i} ]=ts2mat(CEEMD_IMF(:,i), (RATIO(1)+RATIO(2)), max(HORIZON_array), feedbackDelays,ispoint);
%     end
%     
%     for i=1:size(CEEMDAN_IMF, 2)
%         [ CEEMDAN_IMF_trn_data{i}, CEEMDAN_IMF_trn_labels{i}, CEEMDAN_IMF_tst_data{i}, CEEMDAN_IMF_tst_labels{i} ]=ts2mat(CEEMDAN_IMF(:,i), (RATIO(1)+RATIO(2)), max(HORIZON_array), feedbackDelays,ispoint);
%     end
%     
    
    % plot IMF
%     myIMFplot(scaled_data(trnIdx), IMF(trnIdx,:));
    
    % statistical Signi test
%     myIMFstatistical_significance_test(IMF(trnIdx,:)', [0.05, 0.95]);
%     myIMFstatistical_significance_test(CEEMD_IMF(trnIdx,:)', [0.05, 0.95]);
%     myIMFstatistical_significance_test(CEEMDAN_IMF(trnIdx,:)', [0.05, 0.95]);
    
    
    %% Persistent
    % preallocation
    persistent_pred=zeros(size(tst_labels,1), max(HORIZON_array));
    
    PERIOD=48;
    for HORIZON=HORIZON_array
        persistent_pred(:,HORIZON)=scaled_data(tstIdx-PERIOD+HORIZON-1);
    end
    per_residue=tst_labels-repmat(persistent_pred(:,end), 1, max(HORIZON_array));
    [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(persistent_pred(:,end), 1, max(HORIZON_array)), per_residue);
    persistent_RMSE([6 12 18 24 36 48])
    persistent_MASE([6 12 18 24 36 48])
    
end   

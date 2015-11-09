%%
addpath('../dataset/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
HORIZON_array=(1:2:9);
ispoint=0;
numLags=96;
% import data
import_data=importdata('countryenergynsw_2011.csv'); % read data
time=datestr(import_data.textdata(2:end,1));

% correct na entry
load_data=na_correction(import_data.data);

% skip outlier detection
% data=load_data;

%     %% detect outliers
%     IQR_factor=1.5;
%     window=30;
%     threshold=10;
%     outlier_idx1=residual_IQR_outlier( load_data, IQR_factor );
%     outlier_idx2=window_mad_outlier( load_data, window, threshold );
%     % remove outliers
%     outlier_idx=intersect(outlier_idx1, outlier_idx2);
%     data = outlier_correction( load_data, outlier_idx );

% scale
[total_scaled_data, max_data, min_data]=scale_data( load_data,1,0,[],[] );
ts=timeseries(total_scaled_data, time);

% prepare for four typical months
starttime={'2011-01-01 00:00', '2011-05-01 00:00', '2011-07-01 00:00', '2011-10-01 00:00'};
endtime={'2011-01-31 23:30', '2011-05-31 23:30', '2011-07-31 23:30', '2011-10-31 23:30'};
for idx=1:4 % four typical months
    % form a TS
    scaled_ts=getsampleusingtime(ts, starttime{idx}, endtime{idx});
    scaled_data=scaled_ts.data;
    
%     % obtain trn, val, tst index
%     datasize=length(scaled_data);
%     int_ratio=(RATIO*datasize);
%     trnIdx=1:int_ratio(1);
%     valIdx=trnIdx(end)+(1:int_ratio(2));
%     tstIdx=valIdx(end)+(1:int_ratio(3));
    
    % use pacf to determine the feedback Delays.
    par_trn=parcorr(scaled_data, numLags);
    % select the top 10% quantile
    feedbackDelays=find(abs(par_trn)>quantile(abs(par_trn), 0.9));
    %     feedbackDelays=find(abs(par_trn)>1.96/sqrt(length(trnIdx)));
    feedbackDelays=feedbackDelays-1;
    feedbackDelays=feedbackDelays(2:end);
    hiddenLayerSize=2*length(feedbackDelays);
    
    RATIO=[17, 7, 7]/31;
    
    [trn_data,trn_labels,tst_data,tst_labels]=ts2mat(scaled_data, (RATIO(1)+RATIO(2)), max(HORIZON_array), feedbackDelays,ispoint);
    
   trnIdx=1:(size(trn_data,1)+max(feedbackDelays));
   tstIdx=trnIdx(end)+(1:size(tst_data,1));
    
%     [ BPNN_hat, tr ] = myNAR( scaled_data, feedbackDelays', hiddenLayerSize, RATIO, HORIZON_array);
%     
%     trnIdx=tr.trainInd+max(feedbackDelays); % offset feedbackDelays
%     valIdx=tr.valInd+max(feedbackDelays); % offset feedbackDelays
%     tstIdx=tr.testInd+max(feedbackDelays); % offset feedbackDelays
%     
%     RATIO=[17*48-max(feedbackDelays), 7*48, 7*48]/(31*48);
    
%     datasize=length(scaled_data);
%     int_ratio=RATIO*datasize;
%     trnIdx=1:int_ratio(1);
%     valIdx=trnIdx(end)+(1:int_ratio(2));
%     tstIdx=valIdx(end)+(1:int_ratio(3));
    
    % EMD IMF series
    %     options.MAXMODES=5;
    IMF=emd(scaled_data);
    IMF=IMF';
    
    % pacf of imfs
    for i=1:size(IMF, 2)
        par=parcorr(IMF(1:round(size(IMF,1)*RATIO(1)),i),numLags);
        tmp_lag=find(abs(par)>1.96/sqrt(round(length(IMF(:,i))*(RATIO(1)+RATIO(2)))))-1;
        imf_lag{i}=setdiff(tmp_lag,0); % remove 0
    end
%     max_imf_lag=max(cell2mat(imf_lag'));
    
    for i=1:size(IMF, 2)
        [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), (RATIO(1)+RATIO(2)), max(HORIZON_array), feedbackDelays,ispoint);
    end
    
    % plot IMF
    myIMFplot(scaled_data(trnIdx), IMF(trnIdx,:));
    
    % statistical Signi test
    myIMFstatistical_significance_test(IMF(trnIdx,:)', [0.05, 0.95]);
    
    %% Persistent
    PERIOD=48;
    for HORIZON=HORIZON_array
        persistent_predict(:,HORIZON)=scaled_data(tstIdx-PERIOD+HORIZON-1);
    end
%     [trn_data,trn_labels,tst_data,tst_labels]=ts2mat(scaled_data, (RATIO(1)+RATIO(2)), 1, feedbackDelays,1);
    %     trn_labels=scaled_data([trnIdx,valIdx]+1);
    %     tst_labels=scaled_data(tstIdx+1);
%     tst_labels_collection{idx}=tst_labels;
    per_residue=tst_labels-repmat(persistent_predict(:,idx), 1, max(HORIZON_array));
    [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(persistent_predict(:,idx), 1, max(HORIZON_array)), per_residue);
    
    % per_residue=tst_labels-repmat(tst_data(:,end), 1, max(HORIZON_array));
    % [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(tst_data(:,end), 1, max(HORIZON_array)), per_residue);
    % persistent_RMSE=persistent_RMSE(HORIZON_array);
    % persistent_sMAPE=persistent_sMAPE(HORIZON_array);
    % persistent_MASE=persistent_MASE(HORIZON_array);
    
    %% BP-NN
    
    for HORIZON=HORIZON_array
        [ bpnn_pred(:,HORIZON) ] = myANN( trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON), hiddenLayerSize, 'trainbfg');
    end
    [ BPNN_RMSE, BPNN_sMAPE, BPNN_MASE ]=myErrorMeasure(tst_labels, bpnn_pred, per_residue);
    
    
%     % use nn toolbox, type `help nntrain' for more info
%     [ BPNN_hat ] = myNAR( scaled_data, feedbackDelays', hiddenLayerSize, RATIO, HORIZON_array);
%     BPNN_predict(:,idx)=BPNN_hat(tstIdx-max(feedbackDelays));
%     [BPNN_RMSE, BPNN_sMAPE, BPNN_MASE]=myErrorMeasure(tst_labels, BPNN_predict(:,idx), per_residue);
%     
    %
    % ANN_RMSE=zeros(1,size(HORIZON_array,2));
    % ANN_sMAPE=zeros(1,size(HORIZON_array,2));
    % ANN_MASE=zeros(1,size(HORIZON_array,2));
    % counter=1;
    % for HORIZON=HORIZON_array
    %     hiddenlayer=2*size(trn_data,2);
    %     [ANN_predict]=myANN(trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON),hiddenlayer); % hidden layer =10
    %     [ANN_RMSE(counter), ANN_sMAPE(counter), ANN_MASE(counter)]=myErrorMeasure(tst_labels(:,HORIZON), ANN_predict, per_residue(:,HORIZON));
    %     counter=counter+1;
    % end
    
    %% RBFNN
    % use newrb function, type `help newrb' for more info
    goal=0.001; % set goal MSE
    spread_range=0.1:0.1:0.9; % set spread
    
    for HORIZON=HORIZON_array
            [ rbfnn_predict(:,HORIZON) ] = ...
                myRBFNN( trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON), goal, spread_range );
    end
        
%     RBFNN_predict(:,idx)=myRBFNN( trn_data, trn_labels, tst_data, tst_labels, goal, spread_range );
    [RBFNN_RMSE, RBFNN_sMAPE, RBFNN_MASE]=myErrorMeasure(tst_labels, rbfnn_predict, per_residue);
    
    % alternatively faster but performance not better
    % RBFNN_pred2=myANN(trn_data, trn_labels, tst_data, tst_labels, 18 , 'trainlm', 'radbas');
    
    %% EMD-BPNN
    for i=1:size(IMF, 2)
        %         [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
        for HORIZON=HORIZON_array
            [ tmp_emd_ann_pred(:,HORIZON) ] = ...
                myANN( IMF_trn_data{i}, IMF_trn_labels{i}(:,HORIZON), IMF_tst_data{i}, IMF_tst_labels{i}(:,HORIZON), hiddenLayerSize );
        end
        %             IMF_CV_accuracy{i}=tmp_IMF_CV_accuracy;
        %             IMF_best_param{i}=tmp_IMF_best_param;
        ANN_IMF_pred{i}=tmp_emd_ann_pred;
    end
    mat_ANNIMF_pred=cell2mat(ANN_IMF_pred);
    for i=1:max(HORIZON_array)
        emdann_pred(:,i)=sum(mat_ANNIMF_pred(:,i:max(HORIZON_array):end),2);
    end
    [ EMD_BPNN_RMSE, EMD_BPNN_sMAPE, EMD_BPNN_MASE ]=myErrorMeasure(tst_labels, emdann_pred, per_residue);
    
    %% EMD-RBFNN
    for i=2:size(IMF, 2)
        %         [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
        for HORIZON=HORIZON_array
            [ tmp_emd_rbfnn_pred(:,HORIZON) ] = ...
                myRBFNN( IMF_trn_data{i}, IMF_trn_labels{i}(:,HORIZON), IMF_tst_data{i}, IMF_tst_labels{i}(:,HORIZON), goal, spread_range );
        end
        %             IMF_CV_accuracy{i}=tmp_IMF_CV_accuracy;
        %             IMF_best_param{i}=tmp_IMF_best_param;
        RBF_IMF_pred{i}=tmp_emd_rbfnn_pred;
    end
    mat_RBFIMF_pred=cell2mat(RBF_IMF_pred);
    for i=1:max(HORIZON_array)
        emdrbf_pred(:,i)=sum(mat_RBFIMF_pred(:,i:max(HORIZON_array):end),2);
    end
    [ EMD_RBFNN_RMSE, EMD_RBFNN_sMAPE, EMD_RBFNN_MASE ]=myErrorMeasure(tst_labels, emdrbf_pred, per_residue);
    
    
    
    
%     for i=1:size(IMF,2)
%         [ IMF_hat(:,i) ] = myNAR( IMF(:,i), feedbackDelays', hiddenLayerSize, RATIO);
%     end
%     IMF_predict=IMF_hat(tstIdx-max(feedbackDelays),:);
%     EMDANN_predict(:,idx)=sum(IMF_predict,2);
%     [EMDANN_RMSE, EMDANN_sMAPE, EMDANN_MASE]=myErrorMeasure(tst_labels, EMDANN_predict(:,idx), per_residue);
%     
    %
    %
    % EMDANN_predict=zeros(length(tst_labels), max(HORIZON_array));
    % for HORIZON=HORIZON_array
    %     ith_predict=zeros(length(tst_labels), size(IMF,1));
    %     for i=1:size(IMF,2)
    %
    %         [ith_trn_data, ith_trn_labels, ith_tst_data, ith_tst_labels] = ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
    %         hiddenlayer=2*size(ith_trn_data,2);
    %         [tmp_predict]=myANN(ith_trn_data, ith_trn_labels(:,HORIZON), ith_tst_data, ith_tst_labels(:,HORIZON),10); % hidden layer =10
    %         ith_predict(:,i)=tmp_predict(length(tmp_predict)-length(ith_tst_labels)+1:end);
    %     end
    %     EMDANN_predict(:,HORIZON)=sum(ith_predict,2);
    % end
    %
    % [EMD_ANN_RMSE, EMD_ANN_sMAPE, EMD_ANN_MASE]=myErrorMeasure(tst_labels, EMDANN_predict, per_residue);
    % EMD_ANN_RMSE=EMD_ANN_RMSE(HORIZON_array);
    % EMD_ANN_sMAPE=EMD_ANN_sMAPE(HORIZON_array);
    % EMD_ANN_MASE=EMD_ANN_MASE(HORIZON_array);
    
    
    
%     %% EMD-A-ANN
%     inputDelays=feedbackDelays';
%     % hiddenLayerSize2=10;
%     [ EMDAANN_hat ] = myNARX( IMF, scaled_data, inputDelays, hiddenLayerSize, RATIO);
%     EMDAANN_predict(:,idx)=EMDAANN_hat(tstIdx-max(inputDelays));
%     [EMD_AANN_RMSE, EMD_AANN_sMAPE, EMD_AANN_MASE]=myErrorMeasure(tst_labels, EMDAANN_predict(:,idx), per_residue);
%     
    
    
    % % for HORIZON=HORIZON_array
    %     ith_predict=zeros(length(tst_labels), size(IMF,1));
    %     for i=1:size(IMF,2)
    %         step=length(F_COUNT);
    %         [imf_trn_data(:,(step*(i-1)+(1:step))), imf_trn_labels, imf_tst_data(:,(step*(i-1)+(1:step))), imf_tst_labels] = ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
    %     end
    %     counter=1;
    %     for HORIZON=HORIZON_array
    %         hiddenlayer=2*size(imf_trn_data,2);
    %         [EMDAANN_predict(:,counter)]=myANN(imf_trn_data, trn_labels(:,HORIZON), imf_tst_data, tst_labels(:,HORIZON),50); % hidden layer =10
    %         counter=counter+1;
    %     end
    %
    % [EMD_AANN_RMSE, EMD_AANN_sMAPE, EMD_AANN_MASE]=myErrorMeasure(tst_labels(:,HORIZON_array), EMDAANN_predict, per_residue(:,HORIZON_array));
    
    %         %% EMD-SVR
    %     c_range=-3:4;
    %     g_range=(1:0.1:2)/max(feedbackDelays);
    %     d_range=2;
    %     e_range=-4:-1;
    %     kernel=2;
    %
    %     for i=1:size(IMF, 2)
    %         [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), RATIO(1)+RATIO(2), 1, feedbackDelays, 1);
    % %         for HORIZON=HORIZON_array
    %             [ model, tmp_IMF_CV_accuracy, tmp_IMF_best_param ] = ...
    %                 grid_search_SVR( IMF_trn_data{i}, IMF_trn_labels{i}, kernel, c_range, g_range, d_range, e_range );
    %             tmp_IMF_pred=svmpredict( IMF_tst_labels{i}, IMF_tst_data{i}, model );
    % %         end
    %         IMF_CV_accuracy{i}=tmp_IMF_CV_accuracy;
    %         IMF_best_param{i}=tmp_IMF_best_param;
    %         IMF_pred{i}=tmp_IMF_pred;
    %     end
    %     mat_IMF_pred=cell2mat(IMF_pred);
    %
    %     for i=1:max(HORIZON_array)
    %         emd_pred(:,i)=sum(mat_IMF_pred(:,i:max(HORIZON_array):end),2);
    %     end
    %     [ EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE ]=myErrorMeasure(tst_labels, emd_pred, per_residue);
    %
    %
    %
    %%
    result{idx}=[persistent_RMSE, persistent_sMAPE, persistent_MASE
        BPNN_RMSE, BPNN_sMAPE, BPNN_MASE
        RBFNN_RMSE, RBFNN_sMAPE, RBFNN_MASE
        EMD_BPNN_RMSE, EMD_BPNN_sMAPE, EMD_BPNN_MASE
        EMD_RBFNN_RMSE, EMD_RBFNN_sMAPE, EMD_RBFNN_MASE];
   pred{idx}=[persistent_predict, bpnn_pred, rbfnn_predict, emdann_pred, emdrbf_pred];
end
save('result_load_20140223.mat', 'result', 'pred');
% save('result_emdann_20140117_4.mat', 'result', 'tst_labels_collection', 'persistent_predict', 'ANN_predict', 'EMDANN_predict', 'EMDAANN_predict');
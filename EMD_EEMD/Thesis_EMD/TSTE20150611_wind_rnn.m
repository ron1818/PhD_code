%%
%addpath('../dataset/NDBC/');
addpath('../../libsvm-3.11/');
addpath('../../EMD_EEMD/');
addpath('../../src/');

%% read data
HORIZON_array=(1:12);
ispoint=0;
numLags=24;
data_name_array={'full_41004h2011.csv','full_44009h2011.csv','full_46077h2011.csv'};
for data_idx=1:3
    % import data
    import_data=importdata(strcat('data/',data_name_array{data_idx})); % read data
    time=datestr(import_data.textdata(2:end,1));
    
    % correct na entry
    load_data=na_correction(import_data.data);
    % scale
    [total_scaled_data, max_data, min_data]=scale_data( load_data(:,2),1,0,[],[] );
    ts=timeseries(total_scaled_data, time);
    
    % prepare for four typical months
    starttime={'2011-01-01 00:50', '2011-04-01 00:50', '2011-07-01 00:50', '2011-11-01 00:50'};
    endtime={'2011-01-31 23:50', '2011-04-30 23:50', '2011-07-31 23:50', '2011-11-30 23:50'};
    
    % preallocate
    result=cell(1,4);
    pred=cell(1,4);
    for idx=1:4 % four typical months
        % form a TS
        scaled_ts=getsampleusingtime(ts, starttime{idx}, endtime{idx});
        scaled_data=scaled_ts.data;
        
        % use pacf to determine the feedback Delays.
        par_trn=pacf(scaled_data', numLags);
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
        

        %% elman network
        elman_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        
        %% Persistent
        % preallocation
        persistent_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        
        for HORIZON=HORIZON_array
            persistent_pred(:,HORIZON)=scaled_data(tstIdx-HORIZON);
        end
        per_residue=tst_labels-repmat(persistent_pred(:,idx), 1, max(HORIZON_array));
        [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(persistent_pred(:,idx), 1, max(HORIZON_array)), per_residue);
        
        %% SVR
        
        % a=IMF_tst_labels{1,1}+IMF_tst_labels{1,2}+IMF_tst_labels{1,3}+IMF_tst_labels{1,4}+IMF_tst_labels{1,5}+IMF_tst_labels{1,6}+IMF_tst_labels{1,7}+IMF_tst_labels{1,8};
        % -t kernel_type : set type of kernel function (default 2)
        %   0 -- linear: u'*v
        %   1 -- polynomial: (gamma*u'*v + coef0)^degree
        %   2 -- radial basis function: exp(-gamma*|u-v|^2)
        %   3 -- sigmoid: tanh(gamma*u'*v + coef0)
        %   4 -- precomputed kernel (kernel values in training_set_file)
        kernel=2; % rbf
        kernel1=0; % linear
        kernel2=1; % poly
        
        c_range=-3:4;
        g_range=(1:0.1:2)/size(trn_data,2);
        d_range=2;
        e_range=-4:-1;
        
        % preallocation
        svr_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        
        for HORIZON=HORIZON_array
            %attr=numLags-lag{idx}+1; % attribute selected
            [ model, CV_accuracy(HORIZON), best_param(HORIZON,:) ] = ...
                grid_search_SVR( trn_data, trn_labels(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
            svr_pred(:,HORIZON)=svmpredict( tst_labels(:,HORIZON), tst_data, model );
        end
        [ SVR_RMSE, SVR_sMAPE, SVR_MASE ]=myErrorMeasure(tst_labels, svr_pred, per_residue);
        
        
        %% EMD-SVR
        % preallocation
        emd_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        IMF_pred=cell(1,size(IMF,2));
        for i=1:size(IMF, 2)
            tmp_IMF_pred=zeros(size(tst_labels,1), max(HORIZON_array));
            
            for HORIZON=HORIZON_array
                [ model, tmp_IMF_CV_accuracy(HORIZON), tmp_IMF_best_param(HORIZON,:) ] = ...
                    grid_search_SVR( IMF_trn_data{i}, IMF_trn_labels{i}(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
                tmp_IMF_pred(:,HORIZON)=svmpredict( IMF_tst_labels{i}(:,HORIZON), IMF_tst_data{i}, model );
            end
            IMF_CV_accuracy{i}=tmp_IMF_CV_accuracy;
            IMF_best_param{i}=tmp_IMF_best_param;
            IMF_pred{i}=tmp_IMF_pred;
        end
        mat_IMF_pred=cell2mat(IMF_pred);
        
        for i=1:max(HORIZON_array)
            emd_pred(:,i)=sum(mat_IMF_pred(:,i:max(HORIZON_array):end),2);
        end
        
        [ EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE ]=myErrorMeasure(tst_labels, emd_pred, per_residue);
        
        %% EEMD-SVR
        % preallocation
        eemd_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        EEMD_IMF_pred=cell(1,size(IMF,2));
        for i=1:size(EEMD_IMF, 2)
            tmp_EEMD_IMF_pred=zeros(size(tst_labels,1), max(HORIZON_array));
            for HORIZON=HORIZON_array
                [ model, tmp_EEMD_IMF_CV_accuracy(HORIZON), tmp_EEMD_IMF_best_param(HORIZON,:) ] = ...
                    grid_search_SVR( EEMD_IMF_trn_data{i}, EEMD_IMF_trn_labels{i}(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
                tmp_EEMD_IMF_pred(:,HORIZON)=svmpredict( EEMD_IMF_tst_labels{i}(:,HORIZON), EEMD_IMF_tst_data{i}, model );
            end
            EEMD_IMF_CV_accuracy{i}=tmp_EEMD_IMF_CV_accuracy;
            EEMD_IMF_best_param{i}=tmp_EEMD_IMF_best_param;
            EEMD_IMF_pred{i}=tmp_EEMD_IMF_pred;
        end
        mat_EEMD_IMF_pred=cell2mat(EEMD_IMF_pred);
        
        for i=1:max(HORIZON_array)
            eemd_pred(:,i)=sum(mat_EEMD_IMF_pred(:,i:max(HORIZON_array):end),2);
        end
        [ EEMD_SVR_RMSE, EEMD_SVR_sMAPE, EEMD_SVR_MASE ]=myErrorMeasure(tst_labels, eemd_pred, per_residue);
        
        
        %% CEEMD-SVR
        % preallocation
        ceemd_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        CEEMD_IMF_pred=cell(1,size(IMF,2));
        for i=1:size(CEEMD_IMF, 2)
            tmp_CEEMD_IMF_pred=zeros(size(tst_labels,1), max(HORIZON_array));
            for HORIZON=HORIZON_array
                [ model, tmp_CEEMD_IMF_CV_accuracy(HORIZON), tmp_CEEMD_IMF_best_param(HORIZON,:) ] = ...
                    grid_search_SVR( CEEMD_IMF_trn_data{i}, CEEMD_IMF_trn_labels{i}(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
                tmp_CEEMD_IMF_pred(:,HORIZON)=svmpredict( CEEMD_IMF_tst_labels{i}(:,HORIZON), CEEMD_IMF_tst_data{i}, model );
            end
            CEEMD_IMF_CV_accuracy{i}=tmp_CEEMD_IMF_CV_accuracy;
            CEEMD_IMF_best_param{i}=tmp_CEEMD_IMF_best_param;
            CEEMD_IMF_pred{i}=tmp_CEEMD_IMF_pred;
        end
        mat_CEEMD_IMF_pred=cell2mat(CEEMD_IMF_pred);
        
        for i=1:max(HORIZON_array)
            ceemd_pred(:,i)=sum(mat_CEEMD_IMF_pred(:,i:max(HORIZON_array):end),2);
        end
        [ CEEMD_SVR_RMSE, CEEMD_SVR_sMAPE, CEEMD_SVR_MASE ]=myErrorMeasure(tst_labels, ceemd_pred, per_residue);
        
        %% CEEMDAN-SVR
        % preallocation
        ceemdan_pred=zeros(size(tst_labels,1), max(HORIZON_array));
        CEEMDAN_IMF_pred=cell(1,size(IMF,2));
        
        for i=1:size(CEEMDAN_IMF, 2)
            % preallocation
            tmp_CEEMDAN_IMF_pred=zeros(size(tst_labels,1), max(HORIZON_array));
            for HORIZON=HORIZON_array
                [ model, tmp_CEEMDAN_IMF_CV_accuracy(HORIZON), tmp_CEEMDAN_IMF_best_param(HORIZON,:) ] = ...
                    grid_search_SVR( CEEMDAN_IMF_trn_data{i}, CEEMDAN_IMF_trn_labels{i}(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
                tmp_CEEMDAN_IMF_pred(:,HORIZON)=svmpredict( CEEMDAN_IMF_tst_labels{i}(:,HORIZON), CEEMDAN_IMF_tst_data{i}, model );
            end
            CEEMDAN_IMF_CV_accuracy{i}=tmp_CEEMDAN_IMF_CV_accuracy;
            CEEMDAN_IMF_best_param{i}=tmp_CEEMDAN_IMF_best_param;
            CEEMDAN_IMF_pred{i}=tmp_CEEMDAN_IMF_pred;
        end
        mat_CEEMDAN_IMF_pred=cell2mat(CEEMDAN_IMF_pred);
        
        for i=1:max(HORIZON_array)
            ceemdan_pred(:,i)=sum(mat_CEEMDAN_IMF_pred(:,i:max(HORIZON_array):end),2);
        end
        [ CEEMDAN_SVR_RMSE, CEEMDAN_SVR_sMAPE, CEEMDAN_SVR_MASE ]=myErrorMeasure(tst_labels, ceemdan_pred, per_residue);
        %%
        result{idx}=[persistent_RMSE, persistent_sMAPE, persistent_MASE
            SVR_RMSE, SVR_sMAPE, SVR_MASE
            EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE
            EEMD_SVR_RMSE, EEMD_SVR_sMAPE, EEMD_SVR_MASE
            CEEMD_SVR_RMSE, CEEMD_SVR_sMAPE, CEEMD_SVR_MASE
            CEEMDAN_SVR_RMSE, CEEMDAN_SVR_sMAPE, CEEMDAN_SVR_MASE];
        delays{idx}=feedbackDelays;
        pred{idx}=[persistent_pred, svr_pred, emd_pred, ceemdan_pred, eemd_pred, ceemd_pred];
    end
    save(strcat('result_wind_',date,data_name_array{data_idx},'.mat'), 'result', 'pred', 'delays');
end
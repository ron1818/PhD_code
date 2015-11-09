clc;
clear;

%%
%addpath('../dataset/NDBC/');
addpath('../../libsvm-3.11/');
addpath('../');
addpath('../../src/');

%% read data
HORIZON_array=[24,48];
ispoint=1;
numLags=96;
data_name_array={'690160_2005.csv','701195_2005.csv','702757_2005.csv'};
for data_idx=1:3
    % import data
    import_data=importdata(strcat('data/',data_name_array{data_idx})); % read data
    time=datestr(import_data.textdata(2:end,1));
    
    % correct na entry
    load_data=na_correction(import_data.data(:,11));
    % negative to zero
    load_data(load_data<0)=0;
    % scale
    [total_scaled_data, max_data, min_data]=scale_data( load_data,1,0,[],[] );
    ts=timeseries(total_scaled_data, time);
    
    % prepare for four typical seasons
    starttime={'2005-01-01 01:00', '2005-04-01 01:00', '2005-07-01 01:00', '2005-10-01 01:00'};
    endtime={'2005-02-28 24:00', '2005-05-31 24:00', '2005-08-31 24:00', '2005-11-30 24:00'};
    
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
        
        RATIO=[0.4,0.3,0.3];
        
        [trn_data,trn_labels,tst_data,tst_labels]=ts2mat(scaled_data, (RATIO(1)+RATIO(2)), HORIZON_array, feedbackDelays,ispoint);
        tst_sunrise=sign(tst_labels);
        
        trnIdx=1:(size(trn_data,1)+max(feedbackDelays));
        tstIdx=trnIdx(end)+(1:size(tst_data,1));
        
        % EMD IMF series
        IMF=emd(scaled_data);
        IMF=IMF';
        
        % EEMD IMF series
        Nstd=0.5;
        NR=100;
        MaxIter=500;
        [ EEMD_IMF ] = myEEMD( scaled_data,Nstd,NR,MaxIter );
        EEMD_IMF=EEMD_IMF';
        
        % CEEMD IMF series
        [ CEEMD_IMF ] = myCEEMD( scaled_data,Nstd,NR,MaxIter );
        CEEMD_IMF=CEEMD_IMF';
        
        % CEEMDAN IMF series
        NR=100;
        Nstd=0.5;
        MaxIter=500;
        [ CEEMDAN_IMF ] = ceemdan( scaled_data,Nstd,NR,MaxIter );
        CEEMDAN_IMF=CEEMDAN_IMF';
        
        for i=1:size(IMF, 2)
            [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), (RATIO(1)+RATIO(2)), HORIZON_array, feedbackDelays,ispoint);
        end
        
        for i=1:size(EEMD_IMF, 2)
            [ EEMD_IMF_trn_data{i}, EEMD_IMF_trn_labels{i}, EEMD_IMF_tst_data{i}, EEMD_IMF_tst_labels{i} ]=ts2mat(EEMD_IMF(:,i), (RATIO(1)+RATIO(2)), HORIZON_array, feedbackDelays,ispoint);
        end
        
        for i=1:size(CEEMD_IMF, 2)
            [ CEEMD_IMF_trn_data{i}, CEEMD_IMF_trn_labels{i}, CEEMD_IMF_tst_data{i}, CEEMD_IMF_tst_labels{i} ]=ts2mat(CEEMD_IMF(:,i), (RATIO(1)+RATIO(2)), HORIZON_array, feedbackDelays,ispoint);
        end
        
        for i=1:size(CEEMDAN_IMF, 2)
            [ CEEMDAN_IMF_trn_data{i}, CEEMDAN_IMF_trn_labels{i}, CEEMDAN_IMF_tst_data{i}, CEEMDAN_IMF_tst_labels{i} ]=ts2mat(CEEMDAN_IMF(:,i), (RATIO(1)+RATIO(2)), HORIZON_array, feedbackDelays,ispoint);
        end
        
        
        %     % plot IMF
        %     myIMFplot(scaled_data(trnIdx), IMF(trnIdx,:));
        
        % statistical Signi test
        %     myIMFstatistical_significance_test(IMF(trnIdx,:)', [0.05, 0.95]);
        %     myIMFstatistical_significance_test(CEEMD_IMF(trnIdx,:)', [0.05, 0.95]);
        %     myIMFstatistical_significance_test(CEEMDAN_IMF(trnIdx,:)', [0.05, 0.95]);
        
        
        %% Persistent
        % preallocation
        persistent_pred=zeros(size(tst_labels,1), length(HORIZON_array));
        
        PERIOD=24;
        for h=1:length(HORIZON_array)
            persistent_pred(:,h)=scaled_data(tstIdx-max(PERIOD,HORIZON_array(h))-1);
        end
        per_residue=tst_labels-persistent_pred;
        [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure_solar(tst_labels, persistent_pred, per_residue, tst_sunrise);
        
        %% SVR
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
        g_range=(1:size(trn_data,2))/size(trn_data,2);
        d_range=2;
        e_range=-4:-1;
        
        % preallocation
        svr_pred=zeros(size(tst_labels,1),length(HORIZON_array));
        
        for h=1:length(HORIZON_array)
            %attr=numLags-lag{idx}+1; % attribute selected
            [ model, CV_accuracy(h), best_param(h,:) ] = ...
                grid_search_SVR( trn_data, trn_labels(:,h), kernel, c_range, g_range, d_range, e_range );
            svr_pred(:,h)=svmpredict( tst_labels(:,h), tst_data, model );
        end
        [ SVR_RMSE, SVR_sMAPE, SVR_MASE ]=myErrorMeasure_solar(tst_labels, svr_pred, per_residue,tst_sunrise);
        
        
        %% EMD-SVR
        % preallocation
        emd_pred=zeros(size(tst_labels,1), length(HORIZON_array));
        IMF_pred=cell(1,size(IMF,2));
        for i=1:size(IMF, 2)
            tmp_IMF_pred=zeros(size(tst_labels,1), length(HORIZON_array));
            
            for h=1:length(HORIZON_array)
                [ model, tmp_IMF_CV_accuracy(h), tmp_IMF_best_param(h,:) ] = ...
                    grid_search_SVR( IMF_trn_data{i}, IMF_trn_labels{i}(:,h), kernel, c_range, g_range, d_range, e_range );
                tmp_IMF_pred(:,h)=svmpredict( IMF_tst_labels{i}(:,h), IMF_tst_data{i}, model );
            end
            IMF_CV_accuracy{i}=tmp_IMF_CV_accuracy;
            IMF_best_param{i}=tmp_IMF_best_param;
            IMF_pred{i}=tmp_IMF_pred;
        end
        mat_IMF_pred=cell2mat(IMF_pred);
        
        for i=1:length(HORIZON_array)
            emd_pred(:,i)=sum(mat_IMF_pred(:,i:length(HORIZON_array):end),2);
        end
        [ EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE ]=myErrorMeasure_solar(tst_labels, emd_pred, per_residue,tst_sunrise);
        
        %% EEMD-SVR
        % preallocation
        eemd_pred=zeros(size(tst_labels,1), length(HORIZON_array));
        EEMD_IMF_pred=cell(1,size(IMF,2));
        for i=1:size(EEMD_IMF, 2)
            tmp_EEMD_IMF_pred=zeros(size(tst_labels,1), length(HORIZON_array));
            for h=1:length(HORIZON_array)
                [ model, tmp_EEMD_IMF_CV_accuracy(h), tmp_EEMD_IMF_best_param(h,:) ] = ...
                    grid_search_SVR( EEMD_IMF_trn_data{i}, EEMD_IMF_trn_labels{i}(:,h), kernel, c_range, g_range, d_range, e_range );
                tmp_EEMD_IMF_pred(:,h)=svmpredict( EEMD_IMF_tst_labels{i}(:,h), EEMD_IMF_tst_data{i}, model );
            end
            EEMD_IMF_CV_accuracy{i}=tmp_EEMD_IMF_CV_accuracy;
            EEMD_IMF_best_param{i}=tmp_EEMD_IMF_best_param;
            EEMD_IMF_pred{i}=tmp_EEMD_IMF_pred;
        end
        mat_EEMD_IMF_pred=cell2mat(EEMD_IMF_pred);
        
        for i=1:length(HORIZON_array)
            eemd_pred(:,i)=sum(mat_EEMD_IMF_pred(:,i:length(HORIZON_array):end),2);
        end
        [ EEMD_SVR_RMSE, EEMD_SVR_sMAPE, EEMD_SVR_MASE ]=myErrorMeasure_solar(tst_labels, eemd_pred, per_residue,tst_sunrise);
        
        
        %% CEEMD-SVR
        % preallocation
        ceemd_pred=zeros(size(tst_labels,1), length(HORIZON_array));
        CEEMD_IMF_pred=cell(1,size(IMF,2));
        for i=1:size(CEEMD_IMF, 2)
            tmp_CEEMD_IMF_pred=zeros(size(tst_labels,1), length(HORIZON_array));
            for h=1:length(HORIZON_array)
                [ model, tmp_CEEMD_IMF_CV_accuracy(h), tmp_CEEMD_IMF_best_param(h,:) ] = ...
                    grid_search_SVR( CEEMD_IMF_trn_data{i}, CEEMD_IMF_trn_labels{i}(:,h), kernel, c_range, g_range, d_range, e_range );
                tmp_CEEMD_IMF_pred(:,h)=svmpredict( CEEMD_IMF_tst_labels{i}(:,h), CEEMD_IMF_tst_data{i}, model );
            end
            CEEMD_IMF_CV_accuracy{i}=tmp_CEEMD_IMF_CV_accuracy;
            CEEMD_IMF_best_param{i}=tmp_CEEMD_IMF_best_param;
            CEEMD_IMF_pred{i}=tmp_CEEMD_IMF_pred;
        end
        mat_CEEMD_IMF_pred=cell2mat(CEEMD_IMF_pred);
        
        for i=1:length(HORIZON_array)
            ceemd_pred(:,i)=sum(mat_CEEMD_IMF_pred(:,i:length(HORIZON_array):end),2);
        end
        [ CEEMD_SVR_RMSE, CEEMD_SVR_sMAPE, CEEMD_SVR_MASE ]=myErrorMeasure_solar(tst_labels, ceemd_pred, per_residue,tst_sunrise);
        
        %% CEEMDAN-SVR
        % preallocation
        ceemdan_pred=zeros(size(tst_labels,1), length(HORIZON_array));
        CEEMDAN_IMF_pred=cell(1,size(IMF,2));
        
        for i=1:size(CEEMDAN_IMF, 2)
            % preallocation
            tmp_CEEMDAN_IMF_pred=zeros(size(tst_labels,1), length(HORIZON_array));
            for h=1:length(HORIZON_array)
                [ model, tmp_CEEMDAN_IMF_CV_accuracy(h), tmp_CEEMDAN_IMF_best_param(h,:) ] = ...
                    grid_search_SVR( CEEMDAN_IMF_trn_data{i}, CEEMDAN_IMF_trn_labels{i}(:,h), kernel, c_range, g_range, d_range, e_range );
                tmp_CEEMDAN_IMF_pred(:,h)=svmpredict( CEEMDAN_IMF_tst_labels{i}(:,h), CEEMDAN_IMF_tst_data{i}, model );
            end
            CEEMDAN_IMF_CV_accuracy{i}=tmp_CEEMDAN_IMF_CV_accuracy;
            CEEMDAN_IMF_best_param{i}=tmp_CEEMDAN_IMF_best_param;
            CEEMDAN_IMF_pred{i}=tmp_CEEMDAN_IMF_pred;
        end
        mat_CEEMDAN_IMF_pred=cell2mat(CEEMDAN_IMF_pred);
        
        for i=1:length(HORIZON_array)
            ceemdan_pred(:,i)=sum(mat_CEEMDAN_IMF_pred(:,i:length(HORIZON_array):end),2);
        end
        [ CEEMDAN_SVR_RMSE, CEEMDAN_SVR_sMAPE, CEEMDAN_SVR_MASE ]=myErrorMeasure_solar(tst_labels, ceemdan_pred, per_residue,tst_sunrise);
        %%
        result{idx}=[persistent_RMSE, persistent_sMAPE, persistent_MASE
            SVR_RMSE, SVR_sMAPE, SVR_MASE
            EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE
            EEMD_SVR_RMSE, EEMD_SVR_sMAPE, EEMD_SVR_MASE
            CEEMD_SVR_RMSE, CEEMD_SVR_sMAPE, CEEMD_SVR_MASE
            CEEMDAN_SVR_RMSE, CEEMDAN_SVR_sMAPE, CEEMDAN_SVR_MASE];
        delays{idx}=feedbackDelays;
        pred{idx}=[persistent_pred, svr_pred, emd_pred, ceemdan_pred, eemd_pred, ceemd_pred];
        target{idx}=tst_labels;
    end
    save(strcat('result_solar_',date,data_name_array{data_idx},'.mat'), 'result', 'pred', 'delays','target');
end
clc;
clear;
%% addpath
addpath('../dataset/NTU_wind/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
week_idx=17:37;
result_date=date;
HORIZON_array=1:5;
F_COUNT=1:12;
RATIO=[0.7, 0.3];

%%
% result_name=['result_', result_date, '_', num2str(HORIZON), '.mat'];


for idx=week_idx
    % write data
    %         fid=fopen(result_name, 'w');
    
    filename = sprintf('marine_dr_week_%d.csv', idx);
    data = importdata(filename);
    data.textdata(1)=[];
    time = datestr(data.textdata);
    ts = timeseries(data.data(:,1), time, 'name', data.textdata{1,1});
    ts.timeinfo.units='10mins';
    ts.timeinfo.format='ddmmmyyyy hh:mm:ss';
    
    %% interpolate NaN
    ts.data=na_correction(ts.data);
    
%     %% detect outliers
%     IQR_factor=1.5;
%     window=30;
%     threshold=10;
%     outlier_idx1=residual_IQR_outlier( ts.data, IQR_factor );
%     outlier_idx2=window_mad_outlier( ts.data, window, threshold );
%     % remove outliers
%     outlier_idx=intersect(outlier_idx1, outlier_idx2);
%     ts.data = outlier_correction( ts.data, outlier_idx );
%     
    %% split training and testing
    % split
    ispoint=0;
    [scaled_data, max_speed, min_speed]=scale_data( ts.data,1,0,[],[] );
    [ trn_data, trn_labels, tst_data, tst_labels ]=ts2mat(scaled_data, RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
    
    
%     trn_end_idx=round(RATIO(1)*length(ts.data));
%     org_trn_speed=ts.data(1:trn_end_idx);
%     org_tst_speed=ts.data(trn_end_idx+1:end);
%     % scale
%     [trn_speed, max_trn_speed, min_trn_speed]=scale_data( org_trn_speed,1,0,[],[] );
%     tst_speed=scale_data( org_tst_speed,1,1e-3,max_trn_speed, min_trn_speed );
%     
%     % TS to vector
%     [trn_data, trn_labels] = ts2mat(trn_speed, max(HORIZON_array), F_COUNT, ispoint);
%     [tst_data, tst_labels] = ts2mat(tst_speed, max(HORIZON_array), F_COUNT, ispoint);
    
options.MAXMODES=5;
IMF=emd(scaled_data, options);
IMF=IMF';

%     % EMD
%     options.MAXMODES=5;
%     IMF= emd([trn_speed; tst_speed], options);
%     
%     % EMD split
%     IMF_trn=IMF(:,1:trn_end_idx)';
%     IMF_tst=IMF(:,trn_end_idx+1:end)';
%     
%     % TS to vector
%     imf_trn_data=cell(size(IMF_trn,2),1);
%     imf_trn_labels=cell(size(IMF_trn,2),1);
%     imf_tst_data=cell(size(IMF_tst,2),1);
%     imf_tst_labels=cell(size(IMF_tst,2),1);
%     for imf_idx=1:size(IMF_trn,2)
%         [imf_trn_data{imf_idx}, imf_trn_labels{imf_idx}] = ts2mat(IMF_trn(:,imf_idx), max(HORIZON_array), F_COUNT, ispoint);
%         [imf_tst_data{imf_idx}, imf_tst_labels{imf_idx}] = ts2mat(IMF_tst(:,imf_idx), max(HORIZON_array), F_COUNT, ispoint);
%     end
    
    %% Persistent
    per_residue=tst_labels-repmat(tst_data(:,end), 1, length(HORIZON_array));
    [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(tst_data(:,end), 1, length(HORIZON_array)), per_residue);
    
    %     %% ARIMA
    %     max_phi=3;
    %     max_theta=3;
    %     max_d=1;
    %     current_d=0;
    %     [ trn_fit, m, phi, theta, d ] = my_arima( trn_labels, max_phi, max_theta, max_d, current_d );
    %     arima_pred=predict(m, tst_labels, HORIZON);
    %     arima_pred=arima_pred{:};
    %     [ARIMA_RMSE, ARIMA_MAPE, ARIMA_MASE]=myErrorMeasure(tst_labels, arima_pred, per_residue);
    
    %% SVR training and testing
    % -t kernel_type : set type of kernel function (default 2)
    % 	0 -- linear: u'*v
    % 	1 -- polynomial: (gamma*u'*v + coef0)^degree
    % 	2 -- radial basis function: exp(-gamma*|u-v|^2)
    % 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
    % 	4 -- precomputed kernel (kernel values in training_set_file)
    kernel=2; % rbf
    kernel2=1; % poly
%     c_range=-2:4;
%     g_range=(1:0.5:2)/F_COUNT;
%     d_range=2:4;
%     e_range=-5:-1;
    c_range=-1:2;
    g_range=(1:0.5:2)/max(F_COUNT);
    d_range=2;
    e_range=-3:-1;
    for HORIZON=HORIZON_array
        [ model, CV_accuracy(HORIZON), best_param(HORIZON,:) ] = ...
            grid_search_SVR( trn_data, trn_labels(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
        pred(:,HORIZON)=svmpredict( tst_labels(:,HORIZON), tst_data, model );
    end
    [ SVR_RMSE, SVR_sMAPE, SVR_MASE ]=myErrorMeasure(tst_labels, pred, per_residue);
%     [ pred_table, fit_table, onestep_model, CV_accuracy, best_param ] = ...
%         svr_progresive( trn_data, trn_labels(:,1), kernel, c_range, g_range, d_range, e_range, tst_data, tst_labels, max(HORIZON_array) );
%     [ SVR_RMSE, SVR_sMAPE, SVR_MASE ]=myErrorMeasure(tst_labels, pred_table, per_residue);
    
    %% EMD-SVR
    for i=1:size(IMF, 2)
        [ IMF_trn_data{i}, IMF_trn_labels{i}, IMF_tst_data{i}, IMF_tst_labels{i} ]=ts2mat(IMF(:,i), RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
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
    
%% EMD-SVR2 hi-RBF lo-poly
    HF_IMF=sum(IMF(:,1:3),2); % high freq
    LF_IMF=sum(IMF(:,4:end),2); % low freq
    [ HF_IMF_trn_data, HF_IMF_trn_labels, HF_IMF_tst_data, HF_IMF_tst_labels ]=ts2mat(HF_IMF, RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
    [ LF_IMF_trn_data, LF_IMF_trn_labels, LF_IMF_tst_data, LF_IMF_tst_labels ]=ts2mat(LF_IMF, RATIO(1), max(HORIZON_array), F_COUNT, ispoint);
    
        for HORIZON=HORIZON_array
            [ model, HF_CV_accuracy(HORIZON), HF_best_param(HORIZON,:) ] = ...
                grid_search_SVR( HF_IMF_trn_data, HF_IMF_trn_labels(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
            HF_IMF_pred(:,HORIZON)=svmpredict( HF_IMF_tst_labels(:,HORIZON), HF_IMF_tst_data, model );
            
            [ model, LF_CV_accuracy(HORIZON), LF_best_param(HORIZON,:) ] = ...
                grid_search_SVR( LF_IMF_trn_data, LF_IMF_trn_labels(:,HORIZON), kernel2, c_range, g_range, d_range, e_range );
            LF_IMF_pred(:,HORIZON)=svmpredict( LF_IMF_tst_labels(:,HORIZON), LF_IMF_tst_data, model );
        end
        EMD2_pred=HF_IMF_pred+LF_IMF_pred;

    [ EMD2_SVR_RMSE, EMD2_SVR_sMAPE, EMD2_SVR_MASE ]=myErrorMeasure(tst_labels, EMD2_pred, per_residue);

    %% AEMD-SVR EMD as attribute
    
        aemd_trn_data=[IMF_trn_data{1}(:,end), IMF_trn_data{2}(:,end), IMF_trn_data{3}(:,end), IMF_trn_data{4}(:,end), IMF_trn_data{5}(:,end), IMF_trn_data{6}(:,end)];
        aemd_tst_data=[IMF_tst_data{1}(:,end), IMF_tst_data{2}(:,end), IMF_tst_data{3}(:,end), IMF_tst_data{4}(:,end), IMF_tst_data{5}(:,end), IMF_tst_data{6}(:,end)];
        for HORIZON=HORIZON_array
            [ model, tmp_aemd_CV_accuracy(HORIZON), tmp_aemd_best_param(HORIZON,:) ] = ...
                grid_search_SVR( aemd_trn_data, trn_labels(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
                aemd_pred(:,HORIZON)=svmpredict( tst_labels(:,HORIZON), aemd_tst_data, model );
        end
        [ AEMD_SVR_RMSE, AEMD_SVR_sMAPE, AEMD_SVR_MASE ]=myErrorMeasure(tst_labels, aemd_pred, per_residue); 
  
    %% Proposed
    aaemd_trn_data=[IMF_trn_data{1}, IMF_trn_data{2}, IMF_trn_data{3}, IMF_trn_data{4}, IMF_trn_data{5}, IMF_trn_data{6}];
    aaemd_tst_data=[IMF_tst_data{1}, IMF_tst_data{2}, IMF_tst_data{3}, IMF_tst_data{4}, IMF_tst_data{5}, IMF_tst_data{6}];
    for HORIZON=HORIZON_array
            [ model, aaemd_CV_accuracy(HORIZON), aaemd_best_param(HORIZON,:) ] = ...
                grid_search_SVR( aaemd_trn_data, trn_labels(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
                aaemd_pred(:,HORIZON)=svmpredict( tst_labels(:,HORIZON), aaemd_tst_data, model );
        end
        [ AAEMD_SVR_RMSE, AAEMD_SVR_sMAPE, AAEMD_SVR_MASE ]=myErrorMeasure(tst_labels, aaemd_pred, per_residue); 
    
        %% EMD-RBFNN discard first
        goal=0.04;
        spread_range=0.1:0.1:0.5;
%         for HORIZON=HORIZON_array
%             [ rbfnn_pred(:,HORIZON) ] = myRBFNN( trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON), goal, spread_range);
%         end
%         [ RBFNN_RMSE, RBFNN_sMAPE, RBFNN_MASE ]=myErrorMeasure(tst_labels, rbfnn_pred, per_residue);
%         
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
        
    %% EMD-BPNN
    hiddenLayerSize=8;
    for HORIZON=HORIZON_array
        [ ann_pred(:,HORIZON) ] = myANN( trn_data, trn_labels(:,HORIZON), tst_data, tst_labels(:,HORIZON), hiddenLayerSize);
    end
    [ ANN_RMSE, ANN_sMAPE, ANN_MASE ]=myErrorMeasure(tst_labels, ann_pred, per_residue);
    
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
    [ EMD_ANN_RMSE, EMD_ANN_sMAPE, EMD_ANN_MASE ]=myErrorMeasure(tst_labels, emdann_pred, per_residue);
        
    %% Save result
    result{idx}=[persistent_RMSE, persistent_sMAPE, persistent_MASE
        SVR_RMSE, SVR_sMAPE, SVR_MASE
        EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE
        EMD2_SVR_RMSE, EMD2_SVR_sMAPE, EMD2_SVR_MASE
        AEMD_SVR_RMSE, AEMD_SVR_sMAPE, AEMD_SVR_MASE
        AAEMD_SVR_RMSE, AAEMD_SVR_sMAPE, AAEMD_SVR_MASE
        EMD_RBFNN_RMSE, EMD_RBFNN_sMAPE, EMD_RBFNN_MASE
        EMD_ANN_RMSE, EMD_ANN_sMAPE, EMD_ANN_MASE];
    SVR_param{idx}=best_param;
    aaemd_param{idx}=aaemd_best_param;
    

end
save('result_20140108.mat', 'result', 'SVR_param', 'aaemd_param');
%% addpath
addpath('../dataset/NTU_wind/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
week_idx=17:37;
result_name=['result_', date, '.txt'];
% write data
% fid=fopen(result_name, 'w');
for idx=week_idx
    filename = sprintf('marine_dr_week_%d.csv', idx);
    data = importdata(filename);
    data.textdata(1)=[];
    time = datestr(data.textdata);
    ts = timeseries(data.data(:,1), time, 'name', data.textdata{1,1});
    ts.timeinfo.units='10mins';
    ts.timeinfo.format='ddmmmyyyy hh:mm:ss';
    
    %% interpolate NaN
    ts.data=na_correction(ts.data);
    
    %% detect outliers
    IQR_factor=1.5;
    window=30;
    threshold=10;
    outlier_idx1=residual_IQR_outlier( ts.data, IQR_factor );
    outlier_idx2=window_mad_outlier( ts.data, window, threshold );
    % remove outliers
    outlier_idx=intersect(outlier_idx1, outlier_idx2);
    ts.data = outlier_correction( ts.data, outlier_idx );
    
    %% split training and testing
    % split
    HORIZON=1;
    F_COUNT=14;
    ispoint=1;
    RATIO=[5*6*24 2*6*24];
    org_trn_speed=ts.data(1:RATIO(1));
    org_tst_speed=ts.data(RATIO(1)+1:end);
    % scale
    [trn_speed, max_trn_speed, min_trn_speed]=scale_data( org_trn_speed,1,1e-3,[],[] );
    tst_speed=scale_data( org_tst_speed,1,1e-3,max_trn_speed, min_trn_speed );
    
    % determine F_COUNT
%     if exist('parcorr', 'file')
%         [PARCOR, lags, pbounds]=parcorr(trn_speed, min(length(trn_speed)-1, 20)); % pacf
%         F_COUNT=max( [1, find(PARCOR>pbounds(1), 1, 'last')-1, find(PARCOR<pbounds(2), 1, 'last')-1] );
%     else
%         [PARCOR, sig, cil, ciu] = pacf(trn_speed', min(length(trn_speed)-1, 20)); % pacf
%         F_COUNT=max( [1, find(PARCOR>cil(1), 1, 'last')-1, find(PARCOR<ciu(1), 1, 'last')-1] );
%     end
    [trn_data, trn_labels] = ts2mat(trn_speed, HORIZON, F_COUNT, ispoint);
    [tst_data, tst_labels] = ts2mat(tst_speed, HORIZON, F_COUNT, ispoint);
    
    %% Persistent
    per_residue=tst_labels-tst_data(:,end);
    [persistent_RMSE, persistent_MAPE, persistent_MASE]=myErrorMeasure(tst_labels, tst_data(:,end), per_residue);
    
%     %% ARIMA
%     max_phi=3;
%     max_theta=3;
%     max_d=1;
%     current_d=0;
%     [ trn_fit, m, phi, theta, d ] = my_arima( trn_labels, max_phi, max_theta, max_d, current_d );
%     arima_pred=predict(m, tst_labels, HORIZON);
%     arima_pred=arima_pred{:};
%     [ARIMA_RMSE, ARIMA_MAPE, ARIMA_MASE]=myErrorMeasure(tst_labels, arima_pred, per_residue);
%     
%     %% Grey-Mode
%     GM_pred=zeros(length(tst_labels),1);
%     for i=1:length(tst_labels)
%         [ GM_pred(i) ] = GM11( tst_data(i,:), 1 );
%     end
%     [GM_RMSE, GM_MAPE, GM_MASE]=myErrorMeasure(tst_labels, GM_pred, per_residue);
%     
    %% SVR training and testing
    % -s svm_type : set type of SVM (default 0)
    % 	0 -- C-SVC
    % 	1 -- nu-SVC
    % 	2 -- one-class SVM
    % 	3 -- epsilon-SVR
    % 	4 -- nu-SVR
    % -t kernel_type : set type of kernel function (default 2)
    % 	0 -- linear: u'*v
    % 	1 -- polynomial: (gamma*u'*v + coef0)^degree
    % 	2 -- radial basis function: exp(-gamma*|u-v|^2)
    % 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
    % 	4 -- precomputed kernel (kernel values in training_set_file)
    % -d degree : set degree in kernel function (default 3)
    % -g gamma : set gamma in kernel function (default 1/num_features)
    % -r coef0 : set coef0 in kernel function (default 0)
    % -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
    % -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
    % -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
    % -m cachesize : set cache memory size in MB (default 100)
    % -e epsilon : set tolerance of termination criterion (default 0.001)
    % -h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
    % -b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
    % -wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
    % -v n: n-fold cross validation mode
    % -q : quiet mode (no outputs)
    kernel=2; % rbf
    kernel2=1; % linear
    c_range=-2:4;
    g_range=(1:0.5:2)/F_COUNT;
    d_range=2:4;
    e_range=-5:-1;
    [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( trn_data, trn_labels, kernel, c_range, g_range, d_range, e_range );
    
    [ SVR_pred, acc, val ]=svmpredict( tst_labels, tst_data, model );
    [ SVR_RMSE, SVR_MAPE, SVR_MASE ]=myErrorMeasure(tst_labels, SVR_pred, per_residue);
    
    %% EMD-SVR training and testing
    % EMD
    IMF= emd([trn_speed; tst_speed]);
    IMF_trn=IMF(:,1:length(trn_speed))';
    IMF_tst=IMF(:,length(trn_speed)+1:end)';
    
    imf_fit=zeros(length(IMF_trn)-HORIZON-F_COUNT+1, size(IMF_trn,2));
    imf_pred=zeros(length(IMF_tst)-HORIZON-F_COUNT+1, size(IMF_trn,2));
    for imf_idx=1:size(IMF_trn,2)
        [imf_trn_data, imf_trn_labels] = ts2mat(IMF_trn(:,imf_idx), HORIZON, F_COUNT, ispoint);
        [imf_tst_data, imf_tst_labels] = ts2mat(IMF_tst(:,imf_idx), HORIZON, F_COUNT, ispoint);
        
        [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( imf_trn_data, imf_trn_labels, kernel, c_range, g_range, d_range, e_range );
        imf_fit(:,imf_idx) = svmpredict(imf_trn_labels, imf_trn_data, model);
        
        [ imf_pred(:,imf_idx), acc, val ]=svmpredict( imf_tst_labels, imf_tst_data, model );
    end
    
    emd_pred=sum(imf_pred,2);
    [ EMD_SVR_RMSE, EMD_SVR_MAPE, EMD_SVR_MASE ]=myErrorMeasure(tst_labels, emd_pred, per_residue);
    
    % EMD-attributes
    aemd_trn_data=IMF_trn(1:end-HORIZON,:);
    aemd_trn_labels=trn_speed(1+HORIZON:end);
    aemd_tst_data=IMF_tst(1:end-HORIZON,:);
    aemd_tst_labels=tst_speed(1+HORIZON:end);
    
    [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( aemd_trn_data, aemd_trn_labels, kernel2, c_range, g_range, d_range, e_range );
    aemd_fit = svmpredict( aemd_trn_labels, aemd_trn_data, model );
    aemd_fit=aemd_fit(F_COUNT:end);
    [ aemd_pred, acc, val ]=svmpredict( aemd_tst_labels, aemd_tst_data, model );
    aemd_pred=aemd_pred(F_COUNT:end);
    [ AEMD_SVR_RMSE, AEMD_SVR_MAPE, AEMD_SVR_MASE ]=myErrorMeasure(tst_labels, aemd_pred, per_residue);
    
    % EMD+AEMD aggregation
    
    [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( [imf_fit, aemd_fit], trn_labels, kernel2, c_range, g_range, d_range, e_range );
    p_pred1=svmpredict(tst_labels, [imf_pred, aemd_pred], model);
    [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( imf_fit, trn_labels, kernel2, c_range, g_range, d_range, e_range );
    p_pred2=svmpredict(tst_labels, imf_pred, model);
    [P1_SVR_RMSE, P1_SVR_MAPE, P1_SVR_MASE]=myErrorMeasure(tst_labels, p_pred1, per_residue);
    [P2_SVR_RMSE, P2_SVR_MAPE, P2_SVR_MASE]=myErrorMeasure(tst_labels, p_pred2, per_residue);

%     %% EEMD-SVR training and testing
%     % EEMD
%     
%     [EIMFs its]= ceemdan([trn_speed; tst_speed], 0.01, 50,100);
%     EIMF_trn=EIMFs(:,1:length(trn_speed))';
%     EIMF_tst=EIMFs(:,length(trn_speed)+1:end)';
%     eimf_fit=zeros(length(IMF_trn)-HORIZON-F_COUNT+1, size(IMF_trn,2));
%     eimf_pred=zeros(length(EIMF_tst)-HORIZON-F_COUNT+1, size(EIMF_trn,2));
%     for eimf_idx=1:size(EIMF_trn,2)
%         [eimf_trn_data, eimf_trn_labels] = ts2mat(EIMF_trn(:,eimf_idx), HORIZON, F_COUNT, ispoint);
%         [eimf_tst_data, eimf_tst_labels] = ts2mat(EIMF_tst(:,eimf_idx), HORIZON, F_COUNT, ispoint);
%         
%         [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( eimf_trn_data, eimf_trn_labels, kernel, c_range, g_range, d_range, e_range );
%         eimf_fit(:,eimf_idx) = svmpredict(eimf_trn_labels, eimf_trn_data, model);
%         
%         [ eimf_pred(:,eimf_idx), acc, val ]=svmpredict( eimf_tst_labels, eimf_tst_data, model );
%     end
%     eemd_pred=sum(eimf_pred,2);
%     [ EEMD_SVR_RMSE, EEMD_SVR_MAPE, EEMD_SVR_MASE ]=myErrorMeasure(tst_labels, eemd_pred, per_residue);
%     
%     % EEMD-attributes
%     aeemd_trn_data=EIMF_trn(1:end-HORIZON,:);
%     aeemd_trn_labels=trn_speed(1+HORIZON:end,:);
%     aeemd_tst_data=EIMF_tst(1:end-HORIZON,:);
%     aeemd_tst_labels=tst_speed(1+HORIZON:end,:);
%     
%     [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( aeemd_trn_data, aeemd_trn_labels, kernel, c_range, g_range, d_range, e_range );
%     aeemd_fit = svmpredict( aeemd_trn_labels, aemd_trn_data, model );
%     aeemd_fit=aeemd_fit(F_COUNT:end);
%     [ aeemd_pred, acc, val ]=svmpredict( aeemd_tst_labels, aeemd_tst_data, model );
%     aeemd_pred=aeemd_pred(F_COUNT:end);
%     [ AEEMD_SVR_RMSE, AEEMD_SVR_MAPE, AEEMD_SVR_MASE ]=myErrorMeasure(tst_labels, aeemd_pred, per_residue);
%     
%     % EEMD+AEEMD aggregation
%     
%     [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( [eimf_fit, aeemd_fit], trn_labels, kernel, c_range, g_range, d_range, e_range );
%     p_pred3=svmpredict(tst_labels, [eimf_pred, aeemd_pred], model);
%     [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( eimf_fit, trn_labels, kernel, c_range, g_range, d_range, e_range );
%     p_pred4=svmpredict(tst_labels, eimf_pred, model);
%     [P3_SVR_RMSE, P3_SVR_MAPE, P3_SVR_MASE]=myErrorMeasure(tst_labels, p_pred3, per_residue);
%     [P4_SVR_RMSE, P4_SVR_MAPE, P4_SVR_MASE]=myErrorMeasure(tst_labels, p_pred4, per_residue);
    
    %% Save result
    result{idx}=[persistent_RMSE, persistent_MAPE, persistent_MASE
%         ARIMA_RMSE, ARIMA_MAPE, ARIMA_MASE
%         GM_RMSE, GM_MAPE, GM_MASE
        SVR_RMSE, SVR_MAPE, SVR_MASE
        EMD_SVR_RMSE, EMD_SVR_MAPE, EMD_SVR_MASE
        AEMD_SVR_RMSE, AEMD_SVR_MAPE, AEMD_SVR_MASE
%         EEMD_SVR_RMSE, EEMD_SVR_MAPE, EEMD_SVR_MASE
%         AEEMD_SVR_RMSE, AEEMD_SVR_MAPE, AEEMD_SVR_MASE
        P2_SVR_RMSE, P2_SVR_MAPE, P2_SVR_MASE
        P1_SVR_RMSE, P1_SVR_MAPE, P1_SVR_MASE];
        
%         P3_SVR_RMSE, P3_SVR_MAPE, P3_SVR_MASE
%         P4_SVR_RMSE, P4_SVR_MAPE, P4_SVR_MASE];
    
%     fprintf(fid, ' & Data & Wk%d & \\\\\n', idx);
%     fprintf(fid, ' & RMSE & MAPE & MASE \\\\\n');
%     fprintf(fid, 'Persistent & %.4f & %.4f & %.4f \\\\\n', persistent_RMSE, persistent_MAPE, persistent_MASE);
%     fprintf(fid, 'ARIMA & %.4f & %.4f & %.4f \\\\\n', ARIMA_RMSE, ARIMA_MAPE, ARIMA_MASE);
%     fprintf(fid, 'Grey Mode & %.4f & %.4f & %.4f \\\\\n', GM_RMSE, GM_MAPE, GM_MASE);
%     fprintf(fid, 'SVR & %.4f & %.4f & %.4f \\\\\n', SVR_RMSE, SVR_MAPE, SVR_MASE);
%     fprintf(fid, 'EMD-SVR & %.4f & %.4f & %.4f \\\\\n', EMD_SVR_RMSE, EMD_SVR_MAPE, EMD_SVR_MASE);
%     fprintf(fid, 'EMD-Attr-SVR & %.4f & %.4f & %.4f \\\\\n', AEMD_SVR_RMSE, AEMD_SVR_MAPE, AEMD_SVR_MASE);
%     fprintf(fid, 'EEMD-SVR & %.4f & %.4f & %.4f \\\\\n', EEMD_SVR_RMSE, EEMD_SVR_MAPE, EEMD_SVR_MASE);
%     fprintf(fid, 'EEMD-Attr-SVR & %.4f & %.4f & %.4f \\\\\n', AEEMD_SVR_RMSE, AEEMD_SVR_MAPE, AEEMD_SVR_MASE);
%     fprintf(fid, 'EMD-Attr-SVR-SVR & %.4f & %.4f & %.4f \\\\\n', P1_SVR_RMSE, P1_SVR_MAPE, P1_SVR_MASE);
%     fprintf(fid, 'EMD-SVR-SVR & %.4f & %.4f & %.4f \\\\\n', P2_SVR_RMSE, P2_SVR_MAPE, P2_SVR_MASE);
%     fprintf(fid, 'EEMD-Attr-SVR-SVR & %.4f & %.4f & %.4f \\\\\n', P3_SVR_RMSE, P3_SVR_MAPE, P3_SVR_MASE);
%     fprintf(fid, 'EEMD-SVR-SVR & %.4f & %.4f & %.4f \\\\\n', P4_SVR_RMSE, P4_SVR_MAPE, P4_SVR_MASE);
end
save('result_1.mat', 'result'); 
% fclose(fid);
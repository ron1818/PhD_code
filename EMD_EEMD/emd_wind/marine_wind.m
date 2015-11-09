%% addpath
addpath('../dataset/NTU_wind/');
addpath('../libsvm-3.11/');
addpath('../EMD_EEMD/');
addpath('../src/');

%% read data
week_idx=17:37;
result_date=date;
HORIZON_array=1:5;
F_COUNT=12;
RATIO=[0.7, 0.3];

%%
% result_name=['result_', result_date, '_', num2str(HORIZON), '.mat'];
load('result_20131118.mat');

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
    ispoint=0;
    
    trn_end_idx=round(RATIO(1)*length(ts.data));
    org_trn_speed=ts.data(1:trn_end_idx);
    org_tst_speed=ts.data(trn_end_idx+1:end);
    % scale
    [trn_speed, max_trn_speed, min_trn_speed]=scale_data( org_trn_speed,1,0,[],[] );
    tst_speed=scale_data( org_tst_speed,1,1e-3,max_trn_speed, min_trn_speed );
    
    % TS to vector
    [trn_data, trn_labels] = ts2mat(trn_speed, max(HORIZON_array), F_COUNT, ispoint);
    [tst_data, tst_labels] = ts2mat(tst_speed, max(HORIZON_array), F_COUNT, ispoint);
    
    % EMD
    options.MAXMODES=5;
    IMF= emd([trn_speed; tst_speed], options);
    
    % EMD split
    IMF_trn=IMF(:,1:trn_end_idx)';
    IMF_tst=IMF(:,trn_end_idx+1:end)';
    
    % TS to vector
    imf_trn_data=cell(size(IMF_trn,2),1);
    imf_trn_labels=cell(size(IMF_trn,2),1);
    imf_tst_data=cell(size(IMF_tst,2),1);
    imf_tst_labels=cell(size(IMF_tst,2),1);
    for imf_idx=1:size(IMF_trn,2)
        [imf_trn_data{imf_idx}, imf_trn_labels{imf_idx}] = ts2mat(IMF_trn(:,imf_idx), max(HORIZON_array), F_COUNT, ispoint);
        [imf_tst_data{imf_idx}, imf_tst_labels{imf_idx}] = ts2mat(IMF_tst(:,imf_idx), max(HORIZON_array), F_COUNT, ispoint);
    end
    
    %% Persistent
    per_residue=tst_labels-repmat(tst_data(:,end), 1, length(HORIZON_array));
    [persistent_RMSE, persistent_sMAPE, persistent_MASE]=myErrorMeasure(tst_labels, repmat(tst_data(:,end), 1, length(HORIZON_array)), per_residue);
    
    %         %% ARIMA
    %         max_phi=8;
    %         max_theta=8;
    %         max_d=0;
    %         current_d=0;
    %         clear arima_pred
    %         [ trn_fit, m, phi, theta, d ] = my_arima( trn_data(:,end), max_phi, max_theta, max_d, current_d );
    %         for HORIZON=HORIZON_array
    %             arima_pred(HORIZON)=predict(m, tst_data(:,end), HORIZON);
    %         end
    %         arima_pred=cell2mat(arima_pred);
    %         [ARIMA_RMSE, ARIMA_sMAPE, ARIMA_MASE]=myErrorMeasure(tst_labels(6:end,:), arima_pred(6:end,:), per_residue(6:end,:));
    %
    %     %% SVR training and testing
    %     % -t kernel_type : set type of kernel function (default 2)
    %     % 	0 -- linear: u'*v
    %     % 	1 -- polynomial: (gamma*u'*v + coef0)^degree
    %     % 	2 -- radial basis function: exp(-gamma*|u-v|^2)
    %     % 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
    %     % 	4 -- precomputed kernel (kernel values in training_set_file)
    kernel=2; % rbf
    kernel2=0; % linear
    c_range=-2:4;
    g_range=(1:0.5:2)/F_COUNT;
    d_range=2:3;
    e_range=-5:-1;
    %     c_range=-1:2;
    %     g_range=(1:0.5:2)/F_COUNT;
    %     d_range=2;
    %     e_range=-3:-1;
    %% persistent
    [ pred_table, fit_table, onestep_model ] = ...
        svr_progresive( trn_data, trn_labels(:,1), kernel, c_range, g_range, d_range, e_range, tst_data, tst_labels, max(HORIZON_array) );
    [ SVR_RMSE, SVR_sMAPE, SVR_MASE ]=myErrorMeasure(tst_labels, pred_table, per_residue);
    
    %     %% EMD-SVR
    %     %         imf_fit=zeros(length(IMF_trn)-max(HORIZON_array)-F_COUNT+1, size(IMF_trn,2));
    %     %         imf_pred=zeros(length(IMF_tst)-max(HORIZON_array)-F_COUNT+1, size(IMF_trn,2));
    %     imf_fit=cell(1,size(IMF_trn,2));
    %     imf_pred=cell(1,size(IMF_trn,2));
    %
    % %     for HORIZON=HORIZON_array
    %         for imf_idx=1:size(IMF_trn,2)
    %             [ imf_pred{imf_idx}, imf_fit{imf_idx}, onestep_model ] = ...
    %                 svr_progresive( imf_trn_data{imf_idx}, imf_trn_labels{imf_idx}(:,1), kernel, c_range, g_range, d_range, e_range, imf_tst_data{imf_idx}, imf_tst_labels{imf_idx}, max(HORIZON_array) );
    % %             [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
    % %                 grid_search_SVR( imf_trn_data{imf_idx}, imf_trn_labels{imf_idx}(:,HORIZON), kernel, c_range, g_range, d_range, e_range );
    % %             imf_fit{HORIZON,imf_idx} = svmpredict(imf_trn_labels{imf_idx}(:,HORIZON), imf_trn_data{imf_idx}, model);
    % %             [ imf_pred{HORIZON,imf_idx}, acc, val ]=svmpredict( imf_tst_labels{imf_idx}(:,HORIZON), imf_tst_data{imf_idx}, model );
    %         end
    % %     end
    %     imf_fit_mat=cell2mat(imf_fit);
    %     imf_pred_mat=cell2mat(imf_pred);
    %
    %     imf_fit_mat=reshape(imf_fit_mat, size(imf_fit_mat,1)*max(HORIZON_array), size(IMF_trn,2));
    %     imf_pred_mat=reshape(imf_pred_mat, size(imf_pred_mat,1)*max(HORIZON_array), size(IMF_trn,2));
    %     emd_pred=sum(imf_pred_mat,2);
    %     emd_pred=reshape(emd_pred, length(emd_pred)/length(HORIZON_array), length(HORIZON_array));
    %
    %     [ EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE ]=myErrorMeasure(tst_labels, emd_pred, per_residue);
    
    %% EMD-SVR-3
    a=reshape(cell2mat(imf_trn_data(1:3)), size(imf_trn_data{1},1),size(imf_trn_data{1},2)*3); % HI
    b=reshape(cell2mat(imf_trn_data(4:6)), size(imf_trn_data{1},1),size(imf_trn_data{1},2)*3); % LO
    
    a1=reshape(cell2mat(imf_trn_labels(1:3)), size(imf_trn_labels{1},1),size(imf_trn_labels{1},2)*3); % HI
    b1=reshape(cell2mat(imf_trn_labels(4:6)), size(imf_trn_labels{1},1),size(imf_trn_labels{1},2)*3); % LO
    
    c=reshape(cell2mat(imf_tst_data(1:3)), size(imf_tst_data{1},1),size(imf_tst_data{1},2)*3); % HI
    d=reshape(cell2mat(imf_tst_data(4:6)), size(imf_tst_data{1},1),size(imf_tst_data{1},2)*3); % LO
    
    c1=reshape(cell2mat(imf_tst_labels(1:3)), size(imf_tst_labels{1},1),size(imf_tst_labels{1},2)*3); % HI
    d1=reshape(cell2mat(imf_tst_labels(4:6)), size(imf_tst_labels{1},1),size(imf_tst_labels{1},2)*3); % LO
    
    for i=1:F_COUNT
        IMF_trn_HI(:,i)=sum(a(:,(i-1)*3+(1:3)),2);
        IMF_trn_LO(:,i)=sum(b(:,(i-1)*3+(1:3)),2);
        
        IMF_tst_HI(:,i)=sum(c(:,(i-1)*3+(1:3)),2);
        IMF_tst_LO(:,i)=sum(d(:,(i-1)*3+(1:3)),2);
        
    end
    
    for i=1:max(HORIZON_array)
        IMF_trn_labels_HI(:,i)=sum(a1(:,(i-1)*3+(1:3)),2);
        IMF_trn_labels_LO(:,i)=sum(b1(:,(i-1)*3+(1:3)),2);
        IMF_tst_labels_HI(:,i)=sum(c1(:,(i-1)*3+(1:3)),2);
        IMF_tst_labels_LO(:,i)=sum(d1(:,(i-1)*3+(1:3)),2);
    end
    
    [ imf_pred_HI, imf_fit_HI, onestep_model ] = ...
        svr_progresive( IMF_trn_HI, IMF_trn_labels_HI(:,1), kernel, c_range, g_range, d_range, e_range, IMF_tst_HI, IMF_tst_labels_HI, max(HORIZON_array) );
    [ imf_pred_LO, imf_fit_LO, onestep_model ] = ...
        svr_progresive( IMF_trn_LO, IMF_trn_labels_LO(:,1), kernel2, c_range, g_range, d_range, e_range, IMF_tst_LO, IMF_tst_labels_LO, max(HORIZON_array) ); % poly
    for i=1:max(HORIZON_array)
        emd_2_pred(:,i)=sum([imf_pred_HI(:,i), imf_pred_LO(:,i)],2);
    end
    [ EMD_2_SVR_RMSE, EMD_2_SVR_sMAPE, EMD_2_SVR_MASE ]=myErrorMeasure(tst_labels, emd_2_pred, per_residue);
    
    
        %% EMD-SVR2, AEMD-SVR, EMD-AEMD
    %     imf_fit_mat2=reshape(imf_fit_mat, size(trn_labels,1), size(IMF_trn, 2)*length(HORIZON_array));
    %     imf_pred_mat2=reshape(imf_pred_mat, size(tst_labels,1), size(IMF_trn, 2)*length(HORIZON_array));
    %     emd2_fit_RBF=zeros(size(trn_labels));
    %     emd2_pred_RBF=zeros(size(tst_labels));
    %     emd2_fit_Li=zeros(size(trn_labels));
    %     emd2_pred_Li=zeros(size(tst_labels));
    %     p1_pred=zeros(size(tst_labels));
    %     p2_pred=zeros(size(tst_labels));
    %
    %     aemd_fit=zeros(size(trn_labels));
    %     aemd_pred=zeros(size(tst_labels));
    %
        %%
        aemd_trn_data=cell2mat(imf_trn_data');
        aemd_tst_data=cell2mat(imf_tst_data');
        aemd1_trn_data=aemd_trn_data(:,F_COUNT:F_COUNT:end);
        aemd1_tst_data=aemd_tst_data(:,F_COUNT:F_COUNT:end);
        for HORIZON=HORIZON_array
    %         emd2_trn_data=imf_fit_mat2(:,HORIZON:length(HORIZON_array):end);
    %         emd2_tst_data=imf_pred_mat2(:,HORIZON:length(HORIZON_array):end);
    
            hth_trn_labels=trn_labels(:,HORIZON);
            hth_tst_labels=tst_labels(:,HORIZON);
    
    %         [ emd2_model_RBF, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
    %             grid_search_SVR( emd2_trn_data, hth_trn_labels, kernel, c_range, 1/size(emd2_trn_data,2), d_range, e_range );
    %         emd2_fit_RBF(:,HORIZON) = svmpredict(hth_trn_labels, emd2_trn_data, emd2_model_RBF);
    %         [ emd2_pred_RBF(:,HORIZON), acc, val ]=svmpredict( hth_tst_labels, emd2_tst_data, emd2_model_RBF );
    %
    %         [ emd2_model_Li, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
    %             grid_search_SVR( emd2_trn_data, hth_trn_labels, kernel2, c_range, 1/size(emd2_trn_data,2), d_range, e_range );
    %         emd2_fit_Li(:,HORIZON) = svmpredict(hth_trn_labels, emd2_trn_data, emd2_model_Li);
    %         [ emd2_pred_Li(:,HORIZON), acc, val ]=svmpredict( hth_tst_labels, emd2_tst_data, emd2_model_Li );
    %
            [ aemd_model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
                grid_search_SVR( aemd_trn_data, hth_trn_labels, kernel, c_range, 1/size(aemd_trn_data,2), d_range, e_range );
            aemd_fit(:,HORIZON) = svmpredict(hth_trn_labels, aemd_trn_data, aemd_model);
            [ aemd_pred(:,HORIZON), acc, val ]=svmpredict( hth_tst_labels, aemd_tst_data, aemd_model );
            
            [ aemd1_model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
                grid_search_SVR( aemd1_trn_data, hth_trn_labels, kernel, c_range, 1/size(aemd1_trn_data,2), d_range, e_range );
            aemd1_fit(:,HORIZON) = svmpredict(hth_trn_labels, aemd1_trn_data, aemd1_model);
            [ aemd1_pred(:,HORIZON), acc, val ]=svmpredict( hth_tst_labels, aemd1_tst_data, aemd1_model );
    
    %
    %         [ aemd_pred, aemd_fit, onestep_model ] = ...
    %                 svr_progresive( aemd_trn_data, trn_labels(:,1), kernel, c_range, g_range, d_range, e_range, aemd_tst_data, tst_labels, max(HORIZON_array) );
    
    
    
    %         p_trn_data=[ aemd_fit(:,HORIZON), emd2_trn_data];
    %         p_tst_data=[ aemd_pred(:,HORIZON), emd2_tst_data];
    %
    %         [ model_RBF, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
    %             grid_search_SVR( p_trn_data, hth_trn_labels, kernel, c_range, 1/size(p_trn_data,2), d_range, e_range );
    %         [ model_Li, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = ...
    %             grid_search_SVR( p_trn_data, hth_trn_labels, kernel2, c_range, 1/size(p_trn_data,2), d_range, e_range );
    %         [ p1_pred(:,HORIZON), acc, val ]=svmpredict( hth_tst_labels, p_tst_data, model_RBF );
    %         [ p2_pred(:,HORIZON), acc, val ]=svmpredict( hth_tst_labels, p_tst_data, model_Li );
        end
    
    %     [ EMD2_RBF_RMSE, EMD2_RBF_sMAPE, EMD2_RBF_MASE ]=myErrorMeasure(tst_labels, emd2_pred_RBF, per_residue);
    %     [ EMD2_Li_RMSE, EMD2_Li_sMAPE, EMD2_Li_MASE ]=myErrorMeasure(tst_labels, emd2_pred_Li, per_residue);
        [ AEMD_SVR_RMSE, AEMD_SVR_sMAPE, AEMD_SVR_MASE ]=myErrorMeasure(tst_labels, aemd_pred, per_residue);
        [ AEMD1_SVR_RMSE, AEMD1_SVR_sMAPE, AEMD1_SVR_MASE ]=myErrorMeasure(tst_labels, aemd1_pred, per_residue);
    %     [ P1_RMSE, P1_sMAPE, P1_MASE ]=myErrorMeasure(tst_labels, p1_pred, per_residue);
    %     [ P2_RMSE, P2_sMAPE, P2_MASE ]=myErrorMeasure(tst_labels, p2_pred, per_residue);
    
    
    %     %% Save result
    %     result{idx}=[persistent_RMSE, persistent_sMAPE, persistent_MASE
    result_emd{idx}=[EMD_2_SVR_RMSE, EMD_2_SVR_sMAPE, EMD_2_SVR_MASE
            SVR_RMSE, SVR_sMAPE, SVR_MASE
    %         EMD_SVR_RMSE, EMD_SVR_sMAPE, EMD_SVR_MASE
    %         EMD2_RBF_RMSE, EMD2_RBF_sMAPE, EMD2_RBF_MASE
    %         EMD2_Li_RMSE, EMD2_Li_sMAPE, EMD2_Li_MASE
            AEMD_SVR_RMSE, AEMD_SVR_sMAPE, AEMD_SVR_MASE
             AEMD1_SVR_RMSE, AEMD1_SVR_sMAPE, AEMD1_SVR_MASE];
    %         P1_RMSE, P1_sMAPE, P1_MASE
    %         P2_RMSE, P2_sMAPE, P2_MASE];
    
end
save(['result_emd2.mat'], 'result_emd', 'result_arima');
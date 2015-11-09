function [ pred_table, fit_table, onestep_model, CV_accuracy, best_param ] = ...
    svr_progresive( onestep_training_data, onestep_training_labels, kernel, c_range, g_range, d_range, e_range, testing_data, testing_labels, HORIZON )
% one step ahead training
[ onestep_model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( onestep_training_data, onestep_training_labels, kernel, c_range, g_range, d_range, e_range );

% testing
pred_table=zeros(length(testing_labels), HORIZON);
fit_table=zeros(length(onestep_training_labels), HORIZON);
old_tst_data=testing_data;
old_trn_data=onestep_training_data;
dummy_tst_labels=zeros(size(testing_labels,1),1);
dummy_trn_labels=zeros(size(onestep_training_labels,1),1);
for i =1:HORIZON
    ith_tst_data=old_tst_data;
    ith_trn_data=old_trn_data;
    ith_pred = svmpredict(dummy_tst_labels, ith_tst_data, onestep_model); % dummy tst_labels
    ith_fit = svmpredict(dummy_trn_labels, ith_trn_data, onestep_model);
    pred_table(:,i)=ith_pred; % save each prediction iteration
    fit_table(:,i)=ith_fit; % save each prediction iteration
    % rearrange
    old_tst_data=[ith_tst_data(:,2:end), ith_pred];
    old_trn_data=[ith_trn_data(:,2:end), ith_fit];
end
end


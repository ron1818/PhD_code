function [predict, KNN_mat, best_k, best_d]=knn_regression( trn_data, trn_labels, tst_data, k_fold, k_range, d_matrix, w_flag )
% Lall and Sharma (1996), Nearest Neighbor Bootstrap, Water Resources
% Research 32 (3) pp. 679-693.

% d_matrix mxn: m is the possible trials, n is logic selections 1 means the
% feature is selected and 0 means the feature is not selected

% calculate w_matrix based on w_flag
w_matrix=distance_weight_calculation(d_matrix, w_flag);

% define k_d grid for GCV grid search
% k_d_grid=zeros(length(k_range), length(d_range));

% training
% loop wrt k and d
for idx_d=1:size(d_matrix,1) % loop through the trials, can be discrete
    d=d_matrix(idx_d,:); % 0 and 1
    d=logical(d); % to logical
    w=w_matrix(idx_d,:); % pre-defined distance weight
    
    % select trn_data
    cv_trn_data=repmat(w,size(trn_data,1),1).*trn_data;
    cv_trn_data=cv_trn_data(:,d);
    
    % create dist mat
    cv_dist_mat=zeros(size(cv_trn_data,1));
    for tstIdx=1:size(cv_trn_data,1)
        for trnIdx=1:size(cv_trn_data,1)
            cv_dist_mat(tstIdx, trnIdx)=sqrt(sum(cv_trn_data(tstIdx,:)-cv_trn_data(trnIdx,:).^2));
        end
    end
    
    % remove diagonal
    % delete self-distance
    cv_dist_mat(cv_dist_mat==0)=inf;
    
    % k fold CV
    cv_mse=zeros(k_fold, max(k_range)); % pre allocation
    cvo=cvpartition(size(trn_data,1),'k',k_fold);
    for cv=1:cvo.NumTestSets
        trIdx=cvo.training(cv);
        teIdx=cvo.test(cv);
        
        % delete testing, testing is unknown for creating dist matrix
        tmp_cv_dist_mat=cv_dist_mat;
        tmp_cv_dist_mat(:,teIdx)=inf;
        
        % sort ascebd
        [cv_KNN_weight, cv_KNN_mat]=sort(tmp_cv_dist_mat, 2, 'ascend');
        
        for idx_k=1:length(k_range)
            k=k_range(idx_k);
            % determine kernel function
            Kernel=(1./(1:k))/sum(1./(1:k));
            % filter in knn
            kth_cv_KNN_mat=cv_KNN_mat(teIdx,1:k);
            knn_trn_labels=trn_labels(kth_cv_KNN_mat);
            % Kerneled prediction
            cv_pred=Kernel*(knn_trn_labels');
            % error
            cv_error=cv_pred'-trn_labels(teIdx);
            %CV
            cv_mse(cv,k)=mean(cv_error.^2);%/(1-1/sum(1./(1:k)))^2;
        end
        clear knn_trn_labels;
    end
    cv_mse=mean(cv_mse);
    cv_mse(cv_mse==0)=inf;
    [row, col]=find(cv_mse==min(cv_mse)); % col  is k
    k_d_grid(col(1),idx_d)=min(cv_mse);
%     clear cv_mse;
end
k_d_grid(k_d_grid==0)=inf;
% find best k and d
[row, col]=find(k_d_grid==min(min(k_d_grid)));
best_k=row(1);
best_d=col(1);

% testing
d=d_matrix(best_d,:); % 0 and 1
    d=logical(d); % to logical
    w=w_matrix(best_d,:); % pre-defined distance weight
    
    % weight trn_data
    weighted_trn_data=repmat(w,size(trn_data,1),1).*trn_data;
    weighted_trn_data=weighted_trn_data(:,d);
    % weight tst_data
    weighted_tst_data=repmat(w,size(tst_data,1),1).*tst_data;
    weighted_tst_data=weighted_tst_data(:,d);
    
    % create dist mat
    dist_mat=zeros(size(weighted_tst_data,1), size(weighted_trn_data,1));
    for tstIdx=1:size(weighted_tst_data,1)
        for trnIdx=1:size(weighted_trn_data,1)
            dist_mat(tstIdx, trnIdx)=sqrt(sum(weighted_tst_data(tstIdx,:)-weighted_trn_data(trnIdx,:).^2));
        end
    end

% create tst vs trn knn mat
[KNN_weight, KNN_mat]=sort(dist_mat, 2, 'ascend');
KNN_mat=KNN_mat(:,1:best_k);

Kernel=(1./(1:best_k))/sum(1./(1:best_k));
% evaluate test data
knn_trn_labels=trn_labels(KNN_mat);
predict=Kernel*(knn_trn_labels');
predict=predict';
end

function w_matrix=distance_weight_calculation(d_matrix, w_flag)
DEMP=0.0625;
weighted_d_matrix=repmat((1:size(d_matrix,2)),size(d_matrix,1),1);
% w_flag: 1 unweighted;2 linear decay;3 exponential decay
switch w_flag
    case 1
        w_matrix=d_matrix;
    case 2
        w_matrix=sign(d_matrix).*weighted_d_matrix./repmat(sum(weighted_d_matrix,2),1,size(d_matrix,2));
    case 3
        w_matrix=sign(d_matrix).*exp(DEMP*(weighted_d_matrix));
    otherwise
        w_matrix=d_matrix;
end

end


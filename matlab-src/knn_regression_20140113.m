function [predict, KNN_mat, best_k, best_d]=knn_regression( trn_ts, tst_ts, k_range, d_range, w_flag, HORIZON )
% Lall and Sharma (1996), Nearest Neighbor Bootstrap, Water Resources
% Research 32 (3) pp. 679-693.

% define k_d grid for GCV grid search
k_d_grid=zeros(length(k_range), length(d_range));

% loop wrt k and d
% d_counter=1;
for idx_d=1:length(d_range)
    d=d_range(idx_d);
    trn_data=zeros( size(trn_ts,1)-d-HORIZON+1, size(trn_ts,2)*d);
    trn_labels=zeros( size(trn_ts,1)-d-HORIZON+1, size(trn_ts,2));
    % create d-dimensional feature vector Dt
    for i=1:size(trn_ts,2)
        [trn_data(:, ((i-1)*d+1:(i*d))), trn_labels(:,i)] = ts2mat(trn_ts(:,i), HORIZON, d, 1); % ispoint=1
    end
    trn_labels=sum(trn_labels,2);
    
    % w: 1 unweighted;2 linear decay;3 exponential decay
    if w_flag==2 %linear decay
        w=(1:d)/sum((1:d)); % linear decay
    elseif w_flag==3
        w=exp(-0.0625*(d:-1:1)); % exp reduction
    else
        w=ones(1,d); % unweighted
    end
    
    w=repmat(w,1,size(trn_ts,2));
    % distance matrix trn vs trn
    cv_dist_mat=zeros(size(trn_data,1));
    for tstIdx=1:size(trn_data,1)
        for trnIdx=1:size(trn_data,1)
            cv_dist_mat(tstIdx, trnIdx)=sqrt(sum(w.*((trn_data(tstIdx,:)-trn_data(trnIdx,:)).^2)));
        end
    end
    % remove diagonal
    % delete self-distance
    % cv_dist_mat=cv_dist_mat+diag(repmat(inf, 1, size(cv_dist_mat,1)));
    cv_dist_mat(cv_dist_mat==0)=inf;
    
%     k_counter=1;
    for idx_k=1:length(k_range)
        k=k_range(idx_k);
        % determine kernel function
        Kernel=(1./(1:k))/sum(1./(1:k));
        % knn
        [cv_KNN_weight, cv_KNN_mat]=sort(cv_dist_mat, 2, 'ascend');
        cv_KNN_mat=cv_KNN_mat(:,1:k);

        % leave-one-out CV
        knn_trn_labels=trn_labels(cv_KNN_mat);
        cv_pred=Kernel*(knn_trn_labels');
        cv_error=cv_pred'-trn_labels;
        %GCV
        k_d_grid(k,d)=sum(cv_error.^2)/length(trn_labels);%/(1-1/sum(1./(1:k)))^2;
%         k_counter=k_counter+1;
    end
%     d_counter=d_counter+1;
end
k_d_grid(k_d_grid==0)=inf;
% find best k and d
[row, col]=find(k_d_grid==min(min(k_d_grid)));
best_k=row(1);
best_d=col(1);

trn_data=zeros( size(trn_ts,1)-best_d-HORIZON+1, size(trn_ts,2)*best_d);
trn_labels=zeros( size(trn_ts,1)-best_d-HORIZON+1, size(trn_ts,2));
tst_data=zeros( size(tst_ts,1)-best_d-HORIZON+1, size(tst_ts,2)*best_d);
tst_labels=zeros( size(tst_ts,1)-best_d-HORIZON+1, size(tst_ts,2));
% ts2mat
for i=1:size(trn_ts,2)
    [trn_data(:, ((i-1)*best_d+1:(i*best_d))), trn_labels(:,i)] = ts2mat(trn_ts(:,i), HORIZON, best_d, 1); % ispoint=1
end
trn_labels=sum(trn_labels,2);

for i=1:size(tst_ts,2)
    [tst_data(:, ((i-1)*best_d+1:(i*best_d))), tst_labels(:,i)] = ts2mat(tst_ts(:,i), HORIZON, best_d, 1); % ispoint=1
end
% tst_labels=sum(tst_labels,2);

% [trn_data, trn_labels] = ts2mat(trn_ts, HORIZON, best_d, 1); % ispoint=1
% [tst_data, tst_labels] = ts2mat(tst_ts, HORIZON, best_d, 1); % ispoint=1

% w: 1 unweighted;2 linear decay;3 exponential decay
if w_flag==2 %linear decay
    w=(1:best_d)/sum((1:best_d)); % linear decay
elseif w_flag==3
    w=exp(-0.0625*(best_d:-1:1)); % exp reduction
else
    w=ones(1,best_d); % unweighted
end
w=repmat(w,1,size(trn_ts,2));

% distance matrix tst vs trn
dist_mat=zeros(size(tst_data,1), size(trn_data,1));
for tstIdx=1:size(tst_data,1)
    for trnIdx=1:size(trn_data,1)
        dist_mat(tstIdx, trnIdx)=sqrt(sum(w.*((tst_data(tstIdx,:)-trn_data(trnIdx,:)).^2)));
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

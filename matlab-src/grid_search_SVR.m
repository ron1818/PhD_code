function [ model, CV_accuracy, best_param, coarse_score_grid, fine_score_grid ] = grid_search_SVR( training_data, training_labels, kernel, c_range, g_range, d_range, e_range )
% This function is to provide a two pass grid search on SVR
% inputs:
% training_data: X(t), X(t-1), X(t-h)
% training_labels: X(t+h)
% kernel: 0 linear, 1 poly, 2 rbf and 3 sigmoid
% c_range: C=2^c_range
% g_range: G=g_range
% d_range: D=d_range
% e_range: E=10^e_range
% outputs:
% model: SVR model
% CV_accuracy: best MSE accuracy
% best_param: best C G D E
% coarse_score_grid, fine_score_grid: score grids

% preallocation
if kernel == 0 % linear
    score_grid_length=length(c_range)*length(e_range);
    param_c=zeros(length(c_range)*length(e_range),1); % parameter, c
    param_e=zeros(length(c_range)*length(e_range),1); % parameter, e
    param_d=zeros(score_grid_length,1); % dummy
    param_g=zeros(score_grid_length,1); % dummy
    coarse_score_grid=zeros(score_grid_length,1); % accuracy, c e
    counter=1;
    for a=1:length(c_range)
        for c=1:length(e_range)
            param_c(counter)=c_range(a);
            param_e(counter)=e_range(c);
            counter=counter+1;
        end
    end
   
elseif kernel == 1 % poly
    score_grid_length=length(c_range)*length(g_range)*length(d_range)*length(e_range);
    param_c=zeros(score_grid_length,1); % parameter, c
    param_g=zeros(score_grid_length,1); % parameter, g
    param_d=zeros(score_grid_length,1); % parameter, d
    param_e=zeros(score_grid_length,1); % parameter, e 
    coarse_score_grid=zeros(score_grid_length,1); % accuracy, c g d e
    counter=1;
    for a=1:length(c_range)
        for b=1:length(g_range)
            for d=1:length(d_range)
                for c=1:length(e_range)
                    param_c(counter)=c_range(a);
                    param_g(counter)=g_range(b);
                    param_d(counter)=d_range(d);
                    param_e(counter)=e_range(c);
                    counter=counter+1;
                end
            end
        end
    end   
else % rbf or sigmoid
    score_grid_length=length(c_range)*length(g_range)*length(e_range);
    param_c=zeros(score_grid_length,1); % parameter, c
    param_g=zeros(score_grid_length,1); % parameter, g
    param_e=zeros(score_grid_length,1); % parameter, e
    param_d=zeros(score_grid_length,1); % dummy
    coarse_score_grid=zeros(score_grid_length,1); % accuracy, c g e
    counter=1;
    for a=1:length(c_range)
        for b=1:length(g_range)
            for c=1:length(e_range)
                param_c(counter)=c_range(a);
                param_g(counter)=g_range(b);
                param_e(counter)=e_range(c);
                counter=counter+1;
            end
        end
    end
end

% training
% coarse grid search, 5 fold
parfor k=1:score_grid_length
    fprintf('coarse grid search loop: %d\n', k);
    c=param_c(k);
    g=param_g(k);
    d=param_d(k);
    e=param_e(k);
    options=sprintf('-s 3 -t %u -c %f -g %f -d %f -p %f -v 5 -q 1', kernel, 2^c, g, d, 10^e);
    coarse_score_grid(k)=svmtrain(training_labels, training_data, options);
end

% extract best accuracy and best parameter
coarse_best_accuracy=min(coarse_score_grid); % best cv accuracy
[row]=find( coarse_score_grid==coarse_best_accuracy );
coarse_best_c=param_c(row(1)); % best C
coarse_best_g=param_g(row(1)); % best G
coarse_best_d=param_d(row(1)); % best D
coarse_best_e=param_e(row(1)); % best E

% fine grid search
% define fine grid
fine_c_range=coarse_best_c+(-1:0.5:1);
fine_g_range=coarse_best_g+(0:0.5:1);
fine_e_range=coarse_best_e+(-1:0.5:1);

% pre-allocation, avoid D
fine_score_grid_length=length(fine_c_range)*length(fine_g_range)*length(fine_e_range);
fine_score_grid=zeros(fine_score_grid_length,1); % accuracy, c g e
counter=1;
for a=1:length(fine_c_range)
    for b=1:length(fine_g_range)
        for c=1:length(fine_e_range)
            param_c(counter)=fine_c_range(a);
            param_g(counter)=fine_g_range(b);
            param_e(counter)=fine_e_range(c);
            counter=counter+1;
        end
    end
end

% fine grid search, 5-fold
parfor k=1:fine_score_grid_length
    fprintf('fine grid search loop: %d\n', k);
    c=param_c(k);
    g=param_g(k);
    d=coarse_best_d;
    e=param_e(k);
    options=sprintf('-s 3 -t %u -c %f -g %f -d %f -p %f -v 5 -q 1', kernel, 2^c, g, d, 10^e);
    fine_score_grid(k)=svmtrain(training_labels, training_data, options);
end

% best accuracy and param
fine_best_accuracy=min(fine_score_grid); % best cv accuracy
[row]=find( fine_score_grid==fine_best_accuracy );
fine_best_c=param_c(row(1)); % best C
fine_best_g=param_g(row(1)); % best G
fine_best_e=param_e(row(1)); % best E

best_param=[2^fine_best_c, fine_best_g, coarse_best_d, 10^fine_best_e];

% training
options=sprintf('-s 3 -t %u -c %f -g %f -d %f -p %f -q 1', kernel, 2^fine_best_c, fine_best_g, coarse_best_d, 10^fine_best_e);
model=svmtrain(training_labels, training_data, options);
CV_accuracy=fine_best_accuracy;
end





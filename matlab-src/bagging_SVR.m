function [ models, OOB, OOB_idx, fitted, test_accuracy, predicted_labels, best_param ]...
    = bagging_SVR( training_data, training_labels, testing_data, testing_labels, kernel, c_range, g_range, d_range, e_range, bag_numbers )
% bootstraping
if license('checkout', 'Statistics_Toolbox')
    [tmp, bags] = bootstrp(bag_numbers, [], training_labels, training_data);
elseif exist('randi', 'builtin')
    bags=randi(length(training_labels), length(training_labels), bag_numbers);
else
    bags=randint(length(training_labels), bag_numbers, [1 length(training_labels)]);
end
clear tmp;

if kernel == 0 % linear
    % preallocation
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
    % preallocation
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
    % preallocation
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
% coarse grid search
parfor k=1:score_grid_length
    fprintf('coarse grid search loop: %d\n', k);
    c=param_c(k);
    g=param_g(k);
    d=param_d(k);
    e=param_e(k);
    options=sprintf('-s 3 -t %u -c %f -g %f -d %f -p %f -v 5 -q 1', kernel, 2^c, g, d, 10^e);
    coarse_score_grid(k)=svmtrain(training_labels, training_data, options);
end

coarse_best_accuracy=min(coarse_score_grid); % best cv accuracy
[row]=find( coarse_score_grid==coarse_best_accuracy );
coarse_best_c=param_c(row(1)); % best C
coarse_best_g=param_g(row(1)); % best G
coarse_best_d=param_d(row(1)); % best D
coarse_best_e=param_e(row(1)); % best E

%fine grid search
fine_c_range=coarse_best_c+(-1:0.5:1);
fine_g_range=coarse_best_g+(0:0.5:1);
fine_e_range=coarse_best_e+(-1:0.5:1);

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
    

parfor k=1:fine_score_grid_length
    fprintf('fine grid search loop: %d\n', k);
    c=param_c(k);
    g=param_g(k);
    d=coarse_best_d;
    e=param_e(k);
    options=sprintf('-s 3 -t %u -c %f -g %f -d %f -p %f -v 5 -q 1', kernel, 2^c, g, d, 10^e);
    fine_score_grid(k)=svmtrain(training_labels, training_data, options);
end

fine_best_accuracy=min(fine_score_grid); % best cv accuracy
[row]=find( fine_score_grid==fine_best_accuracy );
fine_best_c=param_c(row(1)); % best C
fine_best_g=param_g(row(1)); % best G
fine_best_e=param_e(row(1)); % best E

best_param=[2^fine_best_c, fine_best_g, coarse_best_d, 10^fine_best_e];
% bagging training
% preallocation
models=cell(1,bag_numbers);
fitted=zeros(length(training_labels), bag_numbers);
residuals=zeros(length(training_labels), bag_numbers);
OOB=zeros(3, bag_numbers);
OOB_idx=zeros(length(training_labels), bag_numbers);
for i=1:bag_numbers % each bag
    ith_bag=bags(:,i);
    ith_data=training_data(ith_bag,:);
    ith_labels=training_labels(ith_bag);
    ith_outbag=setdiff((1:length(training_labels)),ith_bag);
    ith_outbag_data=training_data(ith_outbag,:);
    OOB_idx(ith_outbag,i)=1;
    ith_outbag_labels=training_labels(ith_outbag);
    
    options=sprintf('-s 3 -t %u -c %f -g %f -d %f -p %f -q 1', kernel, 2^fine_best_c, fine_best_g, coarse_best_d, 10^fine_best_e);
    models{i}=svmtrain(ith_labels, ith_data, options);
    [pred]=svmpredict(ith_outbag_labels, ith_outbag_data, models{i});
    OOB(1,i)=sqrt(mean((pred - ith_outbag_labels).^2)); % RMSE
%     OOB(2,i)=mean(abs(pred - ith_outbag_labels)); % MAE
    OOB(2,i)=mean(abs((pred - ith_outbag_labels)./ith_outbag_labels))*100; % MAPE
    OOB(3,i)=mase( ith_outbag_labels, pred ); % MASE
    
    fitted(:,i)=svmpredict(training_labels, training_data, models{i});
    residuals(:,i)=fitted(:,i)-training_labels;
end


% bagging testing
base_classifier_accuracy=zeros(3, bag_numbers);
predicted_labels=nan(size(testing_data,1),bag_numbers);
for j=1:bag_numbers% for each bags
    [predicted_labels(:,j)]=svmpredict(testing_labels, testing_data, models{j});
    base_classifier_accuracy(1,j)=sqrt(mean((predicted_labels(:,j) - testing_labels).^2)); % RMSE
    base_classifier_accuracy(2,j)=mean(abs((predicted_labels(:,j) - testing_labels)/mean(testing_labels)))*100; %MAPE
    base_classifier_accuracy(3,j)=mase( testing_labels, predicted_labels(:,j) ); % MASE
end

test_accuracy=zeros(3,1);
aggregated_predicted_labels=mean(predicted_labels,2);
test_accuracy(1)=sqrt(mean((aggregated_predicted_labels - testing_labels).^2)); % RMSE
test_accuracy(2)=mean(abs((aggregated_predicted_labels - testing_labels)./testing_labels))*100; % MAPE
test_accuracy(3)=mase( testing_labels, aggregated_predicted_labels ); % MASE

end




function [trn_idx, tst_idx, val_idx]=split_datetime(datetime, split_method, ratio, isval)
% split_method = 'month': by month
% split_method = 'season': by season
IDX=1:length(datetime);
month_array={'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'};
% season_array={'(01|02|03)', '(04|05|06)', '(07|08|09)', '(10|11|12)'};
season_array={'(12|01|02)', '(03|04|05)', '(06|07|08)', '(09|10|11)'};
if strcmp(split_method, 'month')
    month_idx=cell(length(datetime),length(month_array));
    for i=1:length(month_array)
        month_idx(:,i)=regexp(datetime, ['\d\d\d\d-', month_array{i}]);
    end
    ix=cellfun('isempty',month_idx);
    month_idx(ix)={0};
    month_idx=cell2mat(month_idx);
    month_idx=sparse(month_idx);
    
    trn_idx=cell(length(month_array),1);
    tst_idx=cell(length(month_array),1);
    val_idx=cell(length(month_array),1);
    if isval==0 % no validation
        for i=1:length(month_array)
            ith_ts=IDX(logical(month_idx(:,i)));
            cutting_point=floor(ratio(1)*sum(month_idx(:,i)));
            trn_idx{i}=ith_ts(1:cutting_point);
            tst_idx{i}=ith_ts(cutting_point+1:end);
        end
    else
        for i=1:length(month_array)
            ith_ts=speed(logical(month_idx(:,i)));
            cutting_point_1=floor(ratio(1)*sum(month_idx(:,i)));
            cutting_point_2=floor((ratio(1)+ratio(2))*sum(month_idx(:,i)));
            trn_idx{i}=ith_ts(1:cutting_point_1);
            val_idx{i}=ith_ts(cutting_point_1+1:cutting_point_2);
            tst_idx{i}=ith_ts(cutting_point_2+1:end);
        end
    end
    
elseif strcmp(split_method, 'season')
    season_idx=cell(length(datetime),length(season_array));
    for i=1:length(season_array)
        season_idx(:,i)=regexp(datetime, ['\d\d\d\d-', season_array{i}]);
    end
    ix=cellfun('isempty',season_idx);
    season_idx(ix)={0};
    season_idx=cell2mat(season_idx);
    season_idx=sparse(season_idx);
    
    trn_idx=cell(length(season_array),1);
    tst_idx=cell(length(season_array),1);
    val_idx=cell(length(season_array),1);
   
    if isval==0 % no validation
        for i=1:length(season_array)
            ith_ts=IDX(logical(season_idx(:,i)));
            cutting_point=floor(ratio(1)*sum(season_idx(:,i)));
            trn_idx{i}=ith_ts(1:cutting_point);
            tst_idx{i}=ith_ts(cutting_point+1:end);
        end
    else
        for i=1:length(season_array)
            ith_ts=IDX(logical(season_idx(:,i)));
            cutting_point_1=floor(ratio(1)*sum(season_idx(:,i)));
            cutting_point_2=floor((ratio(1)+ratio(2))*sum(season_idx(:,i)));
            trn_idx{i}=ith_ts(1:cutting_point_1);
            val_idx{i}=ith_ts(cutting_point_1+1:cutting_point_2);
            tst_idx{i}=ith_ts(cutting_point_2+1:end);
        end
    end
end
end

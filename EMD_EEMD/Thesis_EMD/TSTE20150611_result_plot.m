%% wind data analysis
data_name_array={'full_41004h2011.csv','full_44009h2011.csv','full_46077h2011.csv'};
rmse=[];
smape=[];
mase=[];
for data_idx=1:3
    load(strcat('result_wind_','13-Jun-2015',data_name_array{data_idx},'.mat'));
    tmp1=zeros(4,6*12);
    tmp2=zeros(4,6*12);
    tmp3=zeros(4,6*12);
    for month=1:4
        for horizon=1:12
            tmp1(month,(horizon-1)*6+(1:6))=result{1,month}(:,horizon)';
            tmp2(month,(horizon-1)*6+(1:6))=result{1,month}(:,12+horizon)';
            tmp3(month,(horizon-1)*6+(1:6))=result{1,month}(:,24+horizon)';
        end
    end
    rmse=[rmse;tmp1];
    smape=[smape;tmp2];
    mase=[mase;tmp3];
    
end

csvwrite('rmse_wind.tex',rmse);
csvwrite('smape_wind.tex',smape);
csvwrite('mase_wind.tex',mase);




%% solar data analysis
data_name_array={'690160_2005.csv','701195_2005.csv','702757_2005.csv'};
rmse=[];
smape=[];
mase=[];
for data_idx=1:3
    load(strcat('result_solar_','24-Jul-2015',data_name_array{data_idx},'.mat'));
    tmp1=zeros(4,6*2);
    tmp2=zeros(4,6*2);
    tmp3=zeros(4,6*2);
    for month=1:4
        for horizon=1:2
            tmp1(month,(horizon-1)*6+(1:6))=result{1,month}(:,horizon)';
            tmp2(month,(horizon-1)*6+(1:6))=result{1,month}(:,2+horizon)';
            tmp3(month,(horizon-1)*6+(1:6))=result{1,month}(:,4+horizon)';
        end
    end
    rmse=[rmse;tmp1];
    smape=[smape;tmp2];
    mase=[mase;tmp3];
    
end

csvwrite('rmse_solar_24.tex',rmse);
csvwrite('smape_solar_24.tex',smape);
csvwrite('mase_solar_24.tex',mase);
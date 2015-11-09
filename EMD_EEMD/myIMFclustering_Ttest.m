function [ IMF_trn_new, IMF_tst_new, cut_point ] = myIMFclustering_Ttest( IMF_trn, IMF_tst )
%Liu Xingjie, Mi Zengqiang, Li Peng, Mei huawei,
%Study on the Multi-step Forecasting for Wind Speed Based on EMD,
%International Conference on Sustainable Power Generation and Supply (SUPERGEN09), 2009, 1-5

% ttest
S=zeros(size(IMF_trn,1), size(IMF_trn,2)-1);
for sidx=1:size(IMF_trn,2)-1
    S(:,sidx)=sum(IMF_trn(:,1:sidx),2);
end
[h,p]= ttest(S);
cut_point=max([1, find(p==min(p))-1]);
IMF_trn_new=zeros(size(IMF_trn,1), 3); % hi, mid, low

for i=1:3
    IMF_trn_new(:,1)=sum(IMF_trn(:,1:cut_point),2);
    IMF_trn_new(:,2)=sum(IMF_trn(:,cut_point+1:end-1),2);
    IMF_trn_new(:,3)=IMF_trn(:,end);
    
    IMF_tst_new(:,1)=sum(IMF_tst(:,1:cut_point),2);
    IMF_tst_new(:,2)=sum(IMF_tst(:,cut_point+1:end-1),2);
    IMF_tst_new(:,3)=IMF_tst(:,end);
end

end


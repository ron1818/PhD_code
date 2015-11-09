function [ outlier_idx ] = window_mad_outlier( x, window, threshold )
% find ts outliers 
% http://stats.stackexchange.com/questions/1142/simple-algorithm-for-online-outlier-detection-of-a-generic-time-series/1153#1153
% inputs:
% x: 1-D vector
% window: rolling window size
% threshold: mad scale factor
% outputs:
% outlier_idx

% rolling window to calculate med+threshold*MAD
tmp=zeros(length(x)-window+1,1);
for i=1:length(x)-window+1
    ith_window=x(i:i+window-1);
    tmp(i)=median(ith_window)+threshold*mad(ith_window,1);
end
z=[tmp(1)*(1:window-1)';tmp]; % use tmp(1) throughout the initial period
outlier_idx=find(x>z);
end


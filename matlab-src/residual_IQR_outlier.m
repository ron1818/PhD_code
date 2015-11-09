function[ outlier_idx ] = residual_IQR_outlier( x, IQR_factor )
% find ts outliers 
% http://stats.stackexchange.com/questions/1142/simple-algorithm-for-online-outlier-detection-of-a-generic-time-series/1153#1153
% inputs:
% x: 1-D vector
% IQR_factor: inter quantile range factor, default 1.5
% outputs:
% outlier_idx

% loess smoothing of x
xfit=smooth(x, 'loess');
% residual of loess 
residue=x-xfit;
IQR=diff(quantile(residue,[0.25 0.75]));
limits=quantile(residue, [0.25 0.75])+IQR_factor*IQR*[-1 1]; % beyond IQR_Factor range
score = residue<limits(1) | residue>limits(2); % residue outside limits range

outlier_idx=find(score>0);
function [ y ] = outlier_correction( x, outlier_idx )
% use simple MA to correct the outliers
y=x;
outlier_idx=setdiff(outlier_idx, [1 2]); % remove the first two idx if present
for i = 1:length(outlier_idx)
    y(outlier_idx(i))=0.5*(y(outlier_idx(i)-1)+y(outlier_idx(i)-2));
end
end


function [ y ] = na_correction( x )
% use simple MA to correct NaN
y=x;
nan_idx=find(isnan(x));
for i = 1:length(nan_idx)
    y(nan_idx(i))=0.5*(y(nan_idx(i)-1)+y(nan_idx(i)-2));
end
end


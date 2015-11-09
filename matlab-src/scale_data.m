%function takes new max (nmax), new min (nmin), data as input
%and scale all features of data into range [nmin,nmax]. In addition, it
%deletes constant features columns.
%
%Function:
%[tdata, max_data, min_data] = scale_data( data,nmax,nmin,dmax,dmin )
%
%Input:
%data should be samples (row) *features (column) format
%nmax is an integer, default +1
%nmin is an integer, default -1
%dmax is an array, predefined maxminum features
%dmax is an array, predefined minimum features
%
%Output: 
%tdata is scaled data (samples by feature matrix). Each col has min -1 and max 1.
%
%Author: Ashish Anand
%Any bug you can email to anand.ashish@pmail.ntu.edu.sg
%
%Modification: Ren Ye
%University: Nanyang Technological University, Singapore
%Date: 12/19/2011
%ChangLog:
%(1) Add comments
%(2) Change scale method

function [tdata, max_data, min_data] = scale_data( data,nmax,nmin,dmax,dmin )
% dimension of data
[M,N]=size(data);
% 1st pass: find max and min data
if isempty(dmin)
    min_data=min(data);
else
    min_data=dmin;
end
if isempty(dmax)
    max_data=max(data);
else
    max_data=dmax;
end
% define scaled max and min
if isempty(nmin)
    nmin=-1*ones(1,N);
end
if isempty(nmax)
    nmax=1*ones(1,N);
end

feat_diff=max_data-min_data; % difference in feature
tdata=data;
for i=1:N
    tdata(tdata>repmat(max_data,M,1),i)=max_data(i);
    tdata(tdata<repmat(min_data,M,1),i)=min_data(i);
end
% 2nd pass: scale data
dataout = tdata - repmat(min_data,size(tdata,1),1); % offset original data
zero_feat=feat_diff==0; % constant features
non_zero_feat=feat_diff~=0;
feat_diff(zero_feat)=1;
dataout = (dataout./repmat(feat_diff,size(tdata,1),1))*(nmax-nmin); % scale 
tdata(:,non_zero_feat)=dataout(:,non_zero_feat)+nmin;% offset back
end

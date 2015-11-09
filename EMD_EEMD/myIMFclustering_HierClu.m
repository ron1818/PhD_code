function [ IMF_trn_new, IMF_tst_new, cluster_idx ] = myIMFclustering_HierClu( IMF_trn, IMF_tst, min_mode )
%Kuo, C.-Y.; Wei, S.-K. & Tsai, P.-W., 
%Ensemble empirical mode decomposition with supervised cluster analysis,
%Advances in Adaptive Data Analysis, 2013, 5, 1-19

% hierarchical clustering IMF into m clusters and recombine
% D=pdist(IMF', 'correlation');
% squareform(D);

% find minimum distance use correlation
D=linkage(IMF_trn(:,1:(end-1))', 'single', 'correlation');
% plot dendrogram
dendrogram(D);
% clusters
cluster_idx=cluster(D, 'maxclust', min_mode);

% new IMFs
b=unique(cluster_idx); % total clusters
IMF_trn_new=zeros(size(IMF_trn,1), length(b));
IMF_tst_new=zeros(size(IMF_tst,1), length(b));
for i=b'
    IMF_trn_new(:,i)=sum(IMF_trn(:,cluster_idx==i),2);
    IMF_tst_new(:,i)=sum(IMF_tst(:,cluster_idx==i),2);
end

end


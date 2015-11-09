function [ IMF_trn_new, IMF_tst_new, cluster_idx, ApEn_score ] = myIMFscoring_ApEn( IMF_trn, IMF_tst, dim, rnum, tau, min_mode )
% Zhang and Jun, Chaotic TS EEMD-ApEn ESN, Acta Physica Sinica, 2012, 62.

% calculate standard deviations
sd=std(IMF_trn, 0, 2);

ApEn_score = zeros(size(IMF_trn,2), rnum);

% main calculation and display
% figure
for i = 1:rnum
    r = i*0.02;
    for j=1:size(IMF_trn,2)
    ApEn_score(j,i) = ApEn(dim, r*sd(j), IMF_trn(:,j), tau);
    end

end

% find minimum distance use correlation
D=linkage(ApEn_score, 'single', 'correlation');
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
% r = 0.02*(1:rnum);
% plot(r,result(1,:),'o-',r,result(2,:),'o-',r,result(3,:),'o-')
% axis([0 rnum*0.02 0 1.05*max(result(:))])
% legend('sin','chirp','white noise')
% title(['ApEn, m=' num2str(m) ', \tau=' num2str(tau)],'fontsize',14)
% xlabel r

end


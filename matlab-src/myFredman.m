function [F, ave_rank, is_null_reject]=myFredman(x, alpha)
[N k]=size(x);

rank=tiedrank(x');
ave_rank=mean(rank,2);
chi2=12*N/(k*(k+1))*(sum(ave_rank.^2)-k*(k+1)^2/4);
F=(N-1)*chi2./(N*(k-1)-chi2);

dof=[(k-1),(k-1)*(N-1)];
F_critical=finv(1-alpha, dof(1), dof(2))
is_null_reject=F>F_critical;
end

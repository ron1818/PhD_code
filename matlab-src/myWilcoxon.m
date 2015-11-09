function [z, is_null_reject, r_plus, r_minus]=myWilcoxon(x,y)
N=size(x, 1);

diff=x-y;
rank=tiedrank(abs(diff)')';
sign_rank=sign(diff);

r_plus=sum(rank.*(sign_rank>=0));
r_minus=sum(rank.*(sign_rank<=0));

T=min(r_plus, r_minus);

z=(T-0.25*N*(N+1))/sqrt(1/24*N*(N+1)*(2*N+1));

is_null_reject=(z<-1.96)||(z>1.96);
end
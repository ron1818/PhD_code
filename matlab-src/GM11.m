function [ X_pred ] = GM11( X0, HORIZON )
% Kayacan, Ulutas, Kaynak (2010), Grey system TSP, Expert Sys. App. 37, pp.
% 1784-1789
X0=tocolumn(X0);
k=length(X0);
X1=cumsum(X0);
Z1=(X1(1:end-1)+X1(2:end))/2;
B=[-Z1 ones(k-1,1)];
Y=X0(2:end);
A=B\Y;
a=A(1);
b=A(2);
F = @(K) (X0(1)-b/a)*exp(-a*(K))*(1-exp(a));

X_pred=F(k+(1:HORIZON));

end

function y=tocolumn(x)
[m, n]=size(x);
if (m>n)
    y=x;
else
    y=x';
end
end
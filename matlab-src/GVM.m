function [ X_pred ] = GVM( X0, HORIZON )
% Kayacan, Ulutas, Kaynak (2010), Grey system TSP, Expert Sys. App. 37, pp.
% 1784-1789
k=length(X0);
X1=cumsum(X0);
Z1=(X1(1:end-1)+X1(2:end))/2;
B=[-Z1' (Z1.^2)'];
Y=X0(2:end)';
A=B\Y;
a=A(1);
b=A(2);
F = @(K) (a*X0(1)*(a-b*X0(1))./(b*X0(1)+(a-b*X0(1))*exp(a*(K-1)))).*...
    ((1-exp(a))*exp(a*(K-2))./(b*X0(1)+(a-b*X0(1))*exp(a*(K-2))));

X_pred=F(k+(1:HORIZON-1));

end


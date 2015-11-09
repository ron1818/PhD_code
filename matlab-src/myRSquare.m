function [r2]=myRSquare(Y, Y_hat)
r2=1-sum((Y-Y_hat).^2)./sum((Y-repmat(mean(Y),size(Y,1),1)).^2);
end
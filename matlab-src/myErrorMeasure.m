function [ RMSE, sMAPE, MASE, MSE, MAE, MdAE, MdAPE, RMSPE, RMdSPE, MAPE, sMdAPE, MRAE, MdRAE, GMRAE ] = myErrorMeasure( y, y_hat, error_BM )
% Scale dependent: MAE RMSE MAE MdAE
error=y-y_hat;
MSE=mean(error.^2);
RMSE=sqrt(MSE);
MAE=mean(abs(error));
MdAE=median(abs(error));

% Percentage errors: MAPE MdAPE RMSPE RMdSPE
percentage_err=error./y*100;
MAPE=mean(abs(percentage_err));
MdAPE=median(abs(percentage_err));
RMSPE=sqrt(mean(percentage_err.^2));
RMdSPE=sqrt(median(percentage_err.^2));

% Symmetric
sMAPE=mean(200*abs(y-y_hat)./(y+y_hat));
sMdAPE=median(200*abs(y-y_hat)./(y+y_hat));

% Relative: MRAE MdRAE GMRAE
rela_err=error./error_BM;
MRAE=mean(abs(rela_err));
MdRAE=median(abs(rela_err));
GMRAE=exp(mean(log(abs(rela_err))));

% Scaled
n=length(y);
MASE=sum(abs(error))./(n/(n-1)*sum(abs(diff(y))));
end


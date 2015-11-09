function [ RMSE, sMAPE, MASE, MSE, MAE, MdAE, MdAPE, RMSPE, RMdSPE, MAPE, sMdAPE, MRAE, MdRAE, GMRAE ] = myErrorMeasure_solar( y, y_hat, error_BM, sunrise )
% Scale dependent: MAE RMSE MAE MdAE
sunrise(sunrise==0)=NaN;
error=y-y_hat;
nanerror=error.*sunrise;
MSE=nanmean(nanerror.^2);
RMSE=sqrt(MSE);
MAE=nanmean(abs(nanerror));
MdAE=nanmedian(abs(nanerror));

% Percentage errors: MAPE MdAPE RMSPE RMdSPE
percentage_err=error./y*100;
nanpercentage_err=percentage_err.*sunrise;
MAPE=nanmean(abs(nanpercentage_err));
MdAPE=nanmedian(abs(nanpercentage_err));
RMSPE=sqrt(nanmean(nanpercentage_err.^2));
RMdSPE=sqrt(nanmedian(nanpercentage_err.^2));

% Symmetric
sMAPE=nanmean(200*abs((y-y_hat).*sunrise)./((y+y_hat).*sunrise));
sMdAPE=nanmedian(200*abs((y-y_hat).*sunrise)./((y+y_hat).*sunrise));

% Relative: MRAE MdRAE GMRAE
rela_err=error./error_BM;
nanrela_err=rela_err.*sunrise;
MRAE=nanmean(abs(nanrela_err));
MdRAE=nanmedian(abs(nanrela_err));
GMRAE=exp(nanmean(log(abs(nanrela_err))));

% Scaled
n=length(y);
MASE=sum(abs(error))./(n/(n-1)*sum(abs(diff(y))));
end


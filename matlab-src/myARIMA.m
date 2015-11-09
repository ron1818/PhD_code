function [ m, phi, theta, d ] = myARIMA( y, max_phi, max_theta, max_d, current_d, seasonal )

% phi for AR
% theta for MA
% d for delay
% seasonal for seasonality
    
% augmented Dickey-Fuller test URT
h=adftest(y);
if(h || current_d>=max_d ) % do not contain UR, no need difference, or max allowed differencing exceeded
    d=current_d; % return d
    [acf, lags, bounds]=autocorr(y, min(length(y)-1, 10)); % acf
    [pacf, lags, pbounds]=parcorr(y, min(length(y)-1, 10)); % pacf
    max_theta=min( [max_theta, find(acf>bounds(1), 1, 'last')-1, find(acf<bounds(2), 1, 'last')-1] );
    max_phi=min( [max_phi, find(pacf>pbounds(1), 1, 'last')-1, find(pacf<pbounds(2), 1, 'last')-1] );
    max_theta=max([max_theta, 1]);
    max_phi=max([max_phi, 1]);
    
%     z=iddata(y);
    aic_grid=inf(max_phi+1, max_theta+1); % consider 0
    bic_grid=inf(max_phi+1, max_theta+1); % consider 0
    model_cell=cell(max_phi+1, max_theta+1); % consider 0
    for phi=0:max_phi
        for theta=0:max_theta
            if phi==0 && theta==0
                continue;
            else
                % specify the model
                if ~exist( 'seasonal', 'var') % non seasonal
                    if phi==0 % no AR
                        mdl=arima('MALags', theta );
                    elseif theta==0 % no MA
                        mdl=arima('ARLags', phi );
                    else % with AR and MA
                        mdl=arima('MALags', theta, 'ARLags', phi );
                    end
                else % with seasonal
                    if phi==0 % no AR
                        mdl=arima('MALags', theta, 'Seasonality', seasonal, 'SMALags', seasonal );
                    elseif theta==0 % no MA
                        mdl=arima('ARLags', phi, 'Seasonality', seasonal, 'SARLags', seasonal );
                    else % with AR and MA
                        mdl=arima('ARLags', phi, 'MALags', theta, 'Seasonality', seasonal, 'SARLags', seasonal );
                    end
                end
                
                % estimate the model
                [fit, varcov, logL]=estimate(mdl, y);
                model{phi+1, theta+1}=fit;
                [aic_grid(phi+1, theta+1), bic_grid(phi+1, theta+1)]=aicbic(logL,phi+theta,length(y));
            end
        end
    end
    
    % find min bic and corresponding phi, theta
    min_bic=min(min(bic_grid));
    [best_phi, best_theta]=find(bic_grid==min_bic);
    phi=best_phi(1)-1;
    theta=best_theta(1)-1;
    m=model{phi+1, theta+1};
    
%     y_fit=simulate(m, y);
   
else % contain UR, need difference
    if(current_d<max_d)
        d=current_d+1;
        dy=diff(y);
        if ~exist( 'seasonal', 'var') % non seasonal
            [ m, phi, theta, d ]=myARIMA(dy, max_phi, max_theta, max_d, d);
        else % seasonal
            [ m, phi, theta, d ]=myARIMA(dy, max_phi, max_theta, max_d, d, seasonal);
        end
%         y_fit=cumsum([y(1);y_fit]);
    end
end



% plot(y, '.-k');
% hold on
% plot(y_fit, '.-b');
% hold off


end


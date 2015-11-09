%print result
week_idx=17:37;
results_array={'result_1.mat', 'result_2.mat', 'result_4.mat', 'result_6.mat', 'result_8.mat'};

%% Persistent, SVR, EMD-SVR, EMD-A-SVR, EMD-SVR-SVR, EMD-A-SVR-SVR
fid=fopen('result.csv', 'w');
fid2=fopen('wilcoxon.csv', 'w');
fprintf(fid, 'Persistent, SVR, EMD-SVR, EMD3-SVR, EMD-A-SVR, P1, P2, P3, P4, Persistent, SVR, EMD-SVR, EMD3-SVR, EMD-A-SVR, P1, P2, P3, P4, Persistent, SVR, EMD-SVR, EMD3-SVR, EMD-A-SVR, P1, P2, P3, P4\n'); %Title

tmp=1:9;
counter=1;
for res=results_array
    load(res{:});
    for idx=week_idx
        output_RMSE(idx,tmp)=result{idx}(tmp,1)';%RMSE
        output_MAPE(idx,tmp)=result{idx}(tmp,2)';%MAPE
        output_MASE(idx,tmp)=result{idx}(tmp,3)';%MASE
    end
    
    fprintf(fid, '%s\n', res{:}); % horizon
    fprintf(fid, 'RMSE,,,,,,,,,MAPE,,,,,,,,,MASE,,,,,,,,,\n');
    for idx=week_idx
        fprintf(fid, '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,', output_RMSE(idx,:));
        fprintf(fid, '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,', output_MAPE(idx,:));
        fprintf(fid, '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f \n', output_MASE(idx,:));
    end
    
    % Wilcoxon on SVR vs EMD-SVR
    [z(1,counter), n(1,counter), rp(1,counter), rm(1,counter)]=myWilcoxon(output_RMSE(week_idx,2), output_RMSE(week_idx,3));
    [z(2,counter), n(2,counter), rp(2,counter), rm(2,counter)]=myWilcoxon(output_MAPE(week_idx,2), output_MAPE(week_idx,3));
    [z(3,counter), n(3,counter), rp(3,counter), rm(3,counter)]=myWilcoxon(output_MASE(week_idx,2), output_MASE(week_idx,3));
    
    
    % Wilcoxon on EMD-SVR vs P1
    [z1(1,counter), n1(1,counter), rp1(1,counter), rm1(1,counter)]=myWilcoxon(output_RMSE(week_idx,3), output_RMSE(week_idx,6));
    [z1(2,counter), n1(2,counter), rp1(2,counter), rm1(2,counter)]=myWilcoxon(output_MAPE(week_idx,3), output_MAPE(week_idx,6));
    [z1(3,counter), n1(3,counter), rp1(3,counter), rm1(3,counter)]=myWilcoxon(output_MASE(week_idx,3), output_MASE(week_idx,6));
    
    % Wilcoxon on EMD-SVR vs P2
    [z2(1,counter), n2(1,counter), rp2(1,counter), rm2(1,counter)]=myWilcoxon(output_RMSE(week_idx,3), output_RMSE(week_idx,7));
    [z2(2,counter), n2(2,counter), rp2(2,counter), rm2(2,counter)]=myWilcoxon(output_MAPE(week_idx,3), output_MAPE(week_idx,7));
    [z2(3,counter), n2(3,counter), rp2(3,counter), rm2(3,counter)]=myWilcoxon(output_MASE(week_idx,3), output_MASE(week_idx,7));
    
    % Wilcoxon on EMD-SVR vs P3
    [z3(1,counter), n3(1,counter), rp3(1,counter), rm3(1,counter)]=myWilcoxon(output_RMSE(week_idx,3), output_RMSE(week_idx,8));
    [z3(2,counter), n3(2,counter), rp3(2,counter), rm3(2,counter)]=myWilcoxon(output_MAPE(week_idx,3), output_MAPE(week_idx,8));
    [z3(3,counter), n3(3,counter), rp3(3,counter), rm3(3,counter)]=myWilcoxon(output_MASE(week_idx,3), output_MASE(week_idx,8));
    
    % Wilcoxon on EMD-SVR vs P4
    [z4(1,counter), n4(1,counter), rp4(1,counter), rm4(1,counter)]=myWilcoxon(output_RMSE(week_idx,3), output_RMSE(week_idx,9));
    [z4(2,counter), n4(2,counter), rp4(2,counter), rm4(2,counter)]=myWilcoxon(output_MAPE(week_idx,3), output_MAPE(week_idx,9));
    [z4(3,counter), n4(3,counter), rp4(3,counter), rm4(3,counter)]=myWilcoxon(output_MASE(week_idx,3), output_MASE(week_idx,9));
    

    counter=counter+1;
end

fprintf(fid2, 'SVR vs EMD-SVR,,,,,\n');
fprintf(fid2, ', RMSE,,,, MAPE,,,, MASE,,,,\n');
fprintf(fid2, 'Horizon, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection\n');
for j=1:5
    fprintf(fid2, '%s', num2str(j));
    for i=1:3
        if n(i,j)
            tmp='Y';
        else
            tmp='N';
        end
        fprintf(fid2, ', %.2f, %d, %d, %s', z(i,j), rp(i,j), rm(i,j), tmp); % horizon, z, plus rank, minus rank
    end
    fprintf(fid2, '\n');
end
fprintf(fid2, '\n');

fprintf(fid2, 'EMD-SVR vs P1,,,,,\n');
fprintf(fid2, ', RMSE,,,, MAPE,,,, MASE,,,,\n');
fprintf(fid2, 'Horizon, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection\n');
for j=1:5
    fprintf(fid2, '%s', num2str(j));
    for i=1:3
        if n1(i,j)
            tmp='Y';
        else
            tmp='N';
        end
        fprintf(fid2, ', %.2f, %d, %d, %s', z1(i,j), rp1(i,j), rm1(i,j), tmp); % horizon, z, plus rank, minus rank
    end
    fprintf(fid2, '\n');
end
fprintf(fid2, '\n');

fprintf(fid2, 'EMD-SVR vs P2,,,,,\n');
fprintf(fid2, ', RMSE,,,, MAPE,,,, MASE,,,,\n');
fprintf(fid2, 'Horizon, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection\n');
for j=1:5
    fprintf(fid2, '%s', num2str(j));
    for i=1:3
        if n2(i,j)
            tmp='Y';
        else
            tmp='N';
        end
        fprintf(fid2, ', %.2f, %d, %d, %s', z2(i,j), rp2(i,j), rm2(i,j), tmp); % horizon, z, plus rank, minus rank
    end
    fprintf(fid2, '\n');
end
fprintf(fid2, '\n');

fprintf(fid2, 'EMD-SVR vs P3,,,,,\n');
fprintf(fid2, ', RMSE,,,, MAPE,,,, MASE,,,,\n');
fprintf(fid2, 'Horizon, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection\n');
for j=1:5
    fprintf(fid2, '%s', num2str(j));
    for i=1:3
        if n3(i,j)
            tmp='Y';
        else
            tmp='N';
        end
        fprintf(fid2, ', %.2f, %d, %d, %s', z3(i,j), rp3(i,j), rm3(i,j), tmp); % horizon, z, plus rank, minus rank
    end
    fprintf(fid2, '\n');
end
fprintf(fid2, '\n');

fprintf(fid2, 'EMD-SVR vs P4,,,,,\n');
fprintf(fid2, ', RMSE,,,, MAPE,,,, MASE,,,,\n');
fprintf(fid2, 'Horizon, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection, Z, + Rank, - Rank, H0 Rejection\n');
for j=1:5
    fprintf(fid2, '%s', num2str(j));
    for i=1:3
        if n4(i,j)
            tmp='Y';
        else
            tmp='N';
        end
        fprintf(fid2, ', %.2f, %d, %d, %s', z4(i,j), rp4(i,j), rm4(i,j), tmp); % horizon, z, plus rank, minus rank
    end
    fprintf(fid2, '\n');
end
fclose(fid);
fclose(fid2);



% %% Persistent, ARIMA, SVR, EEMD-SVR, EEMD-A-SVR, EEMD-SVR-SVR, EEMD-A-SVR-SVR
% fid=fopen('result_EEMD.csv', 'w');
% fprintf(fid, 'Persistent, ARIMA, SVR, EEMD-SVR, EEMD-A-SVR, EEMD-SVR-SVR, EEMD-A-SVR-SVR\n'); %Title
%  
%  tmp=[1 2 4 7 8 12 11]; 
% for res=results_array
%     load(res{:});
% for idx=week_idx
%     output_RMSE(idx,1:7)=result{idx}(tmp,1)';%RMSE
%     output_MAPE(idx,1:7)=result{idx}(tmp,2)';%MAPE
%     output_MASE(idx,1:7)=result{idx}(tmp,3)';%MASE
% end
% 
% fprintf(fid, '%s\n', res{:}); % horizon
% fprintf(fid, 'RMSE,,,,,,,MAPE,,,,,,,MASE,,,,,,,\n');
% for idx=week_idx
%     fprintf(fid, '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,', output_RMSE(idx,:));
%     fprintf(fid, '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,', output_MAPE(idx,:));
%     fprintf(fid, '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', output_MASE(idx,:));
% end
% end
% fclose(fid);


% %% ALL: Persistent, ARIMA, SVR, EMD-SVR, EEMD-SVR, EMD-A-SVR, EEMD-A-SVR, Proposed(EMD), Proposed(EEMD)
% tmp=[1 2 4 5 7 6 8 11];
% for idx=week_idx
%     output_RMSE(idx,1:8)=result{idx}(tmp,1)';%RMSE
%     output_MAPE(idx,1:8)=result{idx}(tmp,2)';%MAPE
%     output_MASE(idx,1:8)=result{idx}(tmp,3)';%MASE
% end
% 
% %% Wilcoxon on all
% for i=[1 2 3 4 5 6 7 8]
%     for j=[1 2 3 4 5 6 7 8]
%         [z_RMSE(i,j), N_RMSE(i,j), rplus_RMSE(i,j), rminus_RMSE(i,j)]=myWilcoxon(output_RMSE(week_idx,i), output_RMSE(week_idx,j));
%         [z_MAPE(i,j), N_MAPE(i,j), rplus_MAPE(i,j), rminus_MAPE(i,j)]=myWilcoxon(output_MAPE(week_idx,i), output_MAPE(week_idx,j));
%         [z_MASE(i,j), N_MASE(i,j), rplus_MASE(i,j), rminus_MASE(i,j)]=myWilcoxon(output_MASE(week_idx,i), output_MASE(week_idx,j));
%     end
% end
% 
% %% Wilcoxon on SVR, EMD-SVR n EEMD-SVR
% for i=[3 4 5]
%     for j=[3 4 5]
%         [z, n, rp, rm]=myWilcoxon(output_RMSE(week_idx,i), output_RMSE(week_idx,j));
%         
%         if i==j% diagonal: void
%             Mat1(i,j)=NaN;
%         elseif i>j% lower tri: Z
%             Mat1(i,j)= z;
%         else% upper tri: N
%             Mat1(i,j)= rp;
%         end
%     end
% end




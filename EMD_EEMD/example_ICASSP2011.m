% Example of the CEEMDAN performance, used in the work where CEEMDAN was first presented:
% Cite: 
% M.E.TORRES, M.A. COLOMINAS, G. SCHLOTTHAUER, P. FLANDRIN,
%  "A complete Ensemble Empirical Mode decomposition with adaptive noise," 
%  IEEE Int. Conf. on Acoust., Speech and Signal Proc. ICASSP-11, pp. 4144-4147, Prague (CZ)
%
% The code loads the signal ecg.mat
% It is an ECG from the MIT-BIH Normal Sinus Rhythm Database
% available at http://www.physionet.org/cgi-bin/atm/ATM.
% The first 10 seconds of the first channel of the record 16265 has been
% used in the aboved mentioned paper
%
% -------------------------------------------------------------------------
% Date: June 06,2011
% Authors:  Torres ME, Colominas MA, Schlotthauer G, Flandrin P.
% For problems with the code, please contact the authors:  
% To:  macolominas(AT)bioingenieria.edu.ar 
% CC:  metorres(AT)santafe-conicet.gov.ar
% -------------------------------------------------------------------------
%  This version was run on Matlab 7.10.0 (R2010a)
%--------------------------------------------------------------------------

load ('ecg.mat');

Nstd = 0.2;
NR = 500;
MaxIter = 5000;


[modes its]=ceemdan(ecg,0.2,500,5000);
t=1:length(ecg);

[a b]=size(modes);

figure;
subplot(a+1,1,1);
plot(t,ecg);% the ECG signal is in the first row of the subplot
ylabel('ECG')
set(gca,'xtick',[])
axis tight;

for i=2:a
    subplot(a+1,1,i);
    plot(t,modes(i-1,:));
    ylabel (['IMF ' num2str(i-1)]);
    set(gca,'xtick',[])
    xlim([1 length(ecg)])
end;

subplot(a+1,1,a+1)
plot(t,modes(a,:))
ylabel(['IMF ' num2str(a)])
xlim([1 length(ecg)])

figure;
boxplot(its);
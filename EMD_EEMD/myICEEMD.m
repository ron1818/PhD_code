function [ modos its ] = myICEEMD( x,Nstd,NR,MaxIter )
% ICEEMD as in
% J. Lin, "Improved Ensemble Empirical Mode Decomposition and 
% its Applications to Gearbox Fault Signal Processing", International Journal
% of Computer Science Issues, Vol. 9, Issue 6, No 2, November 2012, pp 1694-0814

%¡¡require EMD_EEMD toolbox, myCEEMD function

% decompose with emd (rough
[ IMF ] = emd( x );

% compute average energy of IMFs exclude residue
Ebar=mean(IMF(1:end-1,:).^2,2);

% identify weak sinosindal and weak transcient signals
mean_power=mean_amp.^2;

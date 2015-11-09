clc;
clear;
load demo1.mat
%% define global variable
global C C0 w E_c E_d in_out P_g K N;


%% main function
% Fitness function and numver of variables
fitnessFcn = @my_fun;
constraintFcn = @my_constraint;
fitnessFcn1 = @my_fun1;
constraintFcn1 = @my_constraint1;
numberOfVariables = K*N;

% If decision variables are bounded provide a bound e.g, LB and UB. 
LB=-1*reshape(in_out',1,numberOfVariables);
UB=1*reshape(in_out',1,numberOfVariables);
Bound = [LB;UB]; % If unbounded then Bound = []

x0=reshape(in_out',1,numberOfVariables);
% x_trf = fmincon(fitnessFcn1,x0,[],[],[],[],LB,UB,constraintFcn1,  optimset('Algorithm','trust-region-reflective'));
% x_ip = fmincon(fitnessFcn1,x0,[],[],[],[],LB,UB,constraintFcn1,  optimset('Algorithm','interior-point'));
% x_as = fmincon(fitnessFcn1,x0,[],[],[],[],LB,UB,constraintFcn1,  optimset('Algorithm','active-set'));
% x_sqp = fmincon(fitnessFcn1,x0,[],[],[],[],LB,UB,constraintFcn1,  optimset('Algorithm','sqp'));

% Create an options structure to be passed to GA
% Three options namely 'CreationFcn', 'MutationFcn', and
% 'PopInitRange' are required part of the problem.
options = gaoptimset('CreationFcn',@int_pop,'MutationFcn',@int_mutation, ...
    'PopInitRange',Bound,'Display','diagnose','StallGenL',40,'Generations',150, ...
    'PopulationSize',120,'PlotFcns',{@gaplotbestf,@gaplotbestindiv},'Vectorized','off', 'UseParallel','always');

% [x,fval] = ga(fitnessFcn,numberOfVariables,options);
[x,fval,exitflag] = GA(fitnessFcn,numberOfVariables,[],[],[],[],LB,UB,constraintFcn,options);

val_my_fun(x)
%% plot results
% capacity and initial capacity
figure();
barh([fliplr(C0), fliplr(C-C0)], 'stacked'); xlabel('Capacity (KWh)'); ylabel('EV');
ylim([0,K+1]);xlim([0,max(C)]);

% in_out time
figure();
ha = tight_subplot(K, 1, 0,[.1 .01],[.1 .01]);
for k=1:K
    bar(ha(k), in_out(k,:));
%     ylabel(num2str(k));
    xlim([1 N+1]);
end
set(ha(1:K-1),'XTickLabel',''); set(ha,'YTickLabel','');
xlabel('Time'); ylabel('EV');

% Pg
figure();
bar(P_g);
xlim([1 N+1]); ylim([min(P_g)-1 max(P_g)+1]);
xlabel('Time'); ylabel('Grid Power');


%%
figure();
ha = tight_subplot(K, 1, [0.01 0.01], [.1 .01],[.1 .01]);
delta=reshape(x, N, K)';
for k=1:K
    kth_delta=delta(k,:);
    kth_delta_pos=kth_delta==1;
    kth_delta_neg=-1*(kth_delta==-1);
    bar(ha(k), [kth_delta_pos;kth_delta_neg]', 'stacked');
    yLim(ha(k), [-1 1]);
    xLim(ha(k), [1 N+1]);
    
end
set(ha(1:K-1),'XTickLabel',''); set(ha,'YTickLabel','');
xlabel('Time'); ylabel('EV');

function [ Ebar, Tbar ] = myIMFstatistical_significance_test( IMF, k)
%Z. Wu and N. E. Huang, ``Statistical significance test of IMF'',
%2005, pp 107-127.
%
% Z. Wu and N. E. Huang, ``A Study of the Characteristics of White Noise
% Using the Empirical Mode Decomposition Method'', 2003, pp 1-25.
%
% Input: IMF is the MxN matrix, N is the sample number, M is the IMF count
% Output: Ebar is the average energy and Tbar is the average period by
% counting extrema

% find the number of extrema
for i=1:size(IMF,1)-1 % exclude Residue
    [indmin, indmax, indzer] = extr(IMF(i,:));
    N_extr(i)=length(indmax);
end

Tbar=size(IMF,2)./N_extr; % average period, exclude Residue
Ebar=mean(IMF(1:end-1,:).^2,2); % average energy, exclude Residue

% plot log E vs log T
figure();
% plot theory line ln(E)+ln(T)=0 => y=-x, E=exp(y), T=exp(x)
x=0:0.5:log(max(Tbar))+0.5;
y=-1*x;
loglog(exp(x), exp(y), '-k');

% plot impirical line: ln(E)=0.12-0.934ln(T)=> y=0.12-0.934x, E=exp(y), T=exp(x)
y=0.12-0.934*x;
hold on;
loglog(exp(x), exp(y), '-b');
% plot spread percentiles: ln(E)=-ln(T)+k\sqrt(2/N) exp(ln(T)/2)
% k=-2.236 -0.675 0 0.675 2.236 for 1st 25th 50th 75th 99th percentile
% respectively
% k=[-2.236 -0.675 0.675 2.236]';
k=norminv(k, 0, 1); % percentile convert back by standard norm dist
k=k';
y=-1*repmat(x,length(k),1)+k.*sqrt(2/size(IMF,2))*exp(x./2);
loglog(exp(x), exp(y), '-.');
% plot IMF points
loglog(Tbar, Ebar, '+r');

hold off
xlabel('ln(Tn)');
ylabel('ln(En)');
end

% copied from emd() from EMD_EEMD toolbox
function [indmin, indmax, indzer] = extr(x,t)

if(nargin==1)
    t=1:length(x);
end

m = length(x);

if nargout > 2
    x1=x(1:m-1);
    x2=x(2:m);
    indzer = find(x1.*x2<0);
    
    if any(x == 0)
        iz = find( x==0 );
        indz = [];
        if any(diff(iz)==1)
            zer = x == 0;
            dz = diff([0 zer 0]);
            debz = find(dz == 1);
            finz = find(dz == -1)-1;
            indz = round((debz+finz)/2);
        else
            indz = iz;
        end
        indzer = sort([indzer indz]);
    end
end

d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;


% when two or more successive points have the same value we consider only one extremum in the middle of the constant area
% (only works if the signal is uniformly sampled)

if any(d==0)
    
    imax = [];
    imin = [];
    
    bad = (d==0);
    dd = diff([0 bad 0]);
    debs = find(dd == 1);
    fins = find(dd == -1);
    if debs(1) == 1
        if length(debs) > 1
            debs = debs(2:end);
            fins = fins(2:end);
        else
            debs = [];
            fins = [];
        end
    end
    if length(debs) > 0
        if fins(end) == m
            if length(debs) > 1
                debs = debs(1:(end-1));
                fins = fins(1:(end-1));
                
            else
                debs = [];
                fins = [];
            end
        end
    end
    lc = length(debs);
    if lc > 0
        for k = 1:lc
            if d(debs(k)-1) > 0
                if d(fins(k)) < 0
                    imax = [imax round((fins(k)+debs(k))/2)];
                end
            else
                if d(fins(k)) > 0
                    imin = [imin round((fins(k)+debs(k))/2)];
                end
            end
        end
    end
    
    if length(imax) > 0
        indmax = sort([indmax imax]);
    end
    
    if length(imin) > 0
        indmin = sort([indmin imin]);
    end 
end
end


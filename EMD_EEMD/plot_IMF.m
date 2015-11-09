% tightplot IMF
% require tight_subplot.m
function plot_IMF(IMF)
K=size(IMF,1);
L=size(IMF,2);
imfplot = tight_subplot(K, 1, [0.01 0.01], [.1 .01],[.1 .01]);
for k=1:K-1
    plot(imfplot(k), IMF(k,:));
    xLim(imfplot(k), [0 L]);
    ylabel(imfplot(k), sprintf('IMF$_{%d}$', k), 'Interpreter', 'latex');
    
end
set(imfplot(1:K-1), 'XTickLabel','');
plot(imfplot(K), IMF(K,:));
xLim(imfplot(K), [0 L]);
ylabel(imfplot(K), sprintf('R$_{%d}$', K), 'Interpreter', 'latex');
end
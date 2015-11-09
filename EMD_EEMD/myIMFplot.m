function myIMFplot( X, IMF, trnIdx, tstIdx )
% check if need to plot trn and tst into one graph
if ~exist('trnIdx', 'var')
    % initialize sup plot numbers,
    % require tight_subplot function
    ha = tight_subplot(size(IMF,2)+1, 1, [0.01 0.01], [.1 .01],[.1 .01]);
    
    plot(ha(1), X); % plot original TS
    set(ha(1), 'xtick', []);
    set(ha(1), 'xticklabel', []);
    ylabel(ha(1), 'org');
    for k=1:size(IMF,2)
        plot(ha(k+1), IMF(:,k));
        ylim(ha(k+1), [min(IMF(:,k)), max(IMF(:,k))]);
        if(k<size(IMF,2))
            set(ha(k+1), 'xtick', []);
            set(ha(k+1), 'xticklabel', []);
            ylabel(ha(k+1), sprintf('c%d', k));
        else
            ylabel(ha(k+1), sprintf('r%d', k-1));
        end
    end
else % need to split
    trnX=X(trnIdx);
    trnIMF=IMF(trnIdx,:);
    tstX=X(tstIdx);
    tstIMF=IMF(tstIdx,:);
    
    
    ha = tight_subplot(size(IMF,2)+1, 1, [0.01 0.01], [.1 .01],[.1 .01]);
    
    plot(ha(1), trnIdx, trnX, '-', tstIdx, tstX, '--'); % plot original TS, solid for trn, dashed for tst
    set(ha(1), 'xtick', []);
    set(ha(1), 'xticklabel', []);
    ylabel(ha(1), 'org');
    for k=1:size(IMF,2)
        plot(ha(k+1), trnIdx, trnIMF(:,k), '-', tstIdx, tstIMF(:,k), '--');
        ylim(ha(k+1), [min(IMF(:,k)), max(IMF(:,k))]);
        if(k<size(IMF,2))
            set(ha(k+1), 'xtick', []);
            set(ha(k+1), 'xticklabel', []);
            ylabel(ha(k+1), sprintf('c%d', k));
        else
            ylabel(ha(k+1), sprintf('r%d', (k-1)));
        end
    end
end

end

function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   going row-wise as in
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end

if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh;

ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
end


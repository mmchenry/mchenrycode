function plotFreqResp(varargin)
% Accewpts a variable number of c structures to compare bode plots
% plotFreqResp(c) - simply plots the freq response fro one c
% plotFreqResp(Cs,clr) - plots variable number of c, requires color matrix
% plotFreqResp(Cs,clr,...) - give a string describing each input


Cs       = varargin{1};

if nargin > 1
    clr = varargin{2};
    if nargin > 2
        lins = varargin{3};
    else
        lins = {{'-'}};
    end
else
    clr = [0 0 1];
    lins = {{'-'}};
end

if length(Cs) > 1
    if (length(Cs)~=size(clr,1))
        error('color matrix needs to have one row per freq. resp.')
    end 
end



for ii = 1:length(Cs)
    c = Cs{ii}(1);
    
    if ~isfield(c,'sensitivity')
        error('Need to run calcFreqResponse first');
    end

    subplot(2,1,1)
    
    if length(lins) > 1
        cLine = lins{ii}{1};
    else
        cLine = lins{1}{1};
    end
    
    h = loglog((c.freqs),(c.sensitivity),cLine);
            % log10(c.peak_freq),log10(c.peak_sense),'+');
    set(h,'Color',clr(ii,:))
    hold on

 %   h = plot();
  %  set(h,'MarkerEdgeColor',clr(ii,:))

   
    subplot(2,1,2)

    h = semilogx((c.freqs),c.phase,cLine);
    hold on
    
    set(h,'Color',clr(ii,:))
    
    %axis square
    

end


subplot(2,1,1)
ylabel('Sensitivity (s)')
xlabel('Freq. (Hz)');
%axis equal
%tmp = get(gca,'YTick');
%set(gca,'YTick',[tmp(1):tmp(end)]);
%tmp = get(gca,'XTick');
%set(gca,'XTick',[tmp(1):tmp(end)]);
%grid on
modPlot

subplot(2,1,2)
ylabel('Phase (deg)');
xlabel('Freq. (Hz)');
%grid on
modPlot
    
if nargin > 3
    str{1} = varargin{4};
    if nargin > 4
        for jj = 1:nargin-4
            str{jj+1} = varargin{jj+4};
        end
    end
    subplot(2,1,1)
    legend(str)
end

set(gcf,'Color',[1 1 1])

% set(h,'MarkerSize',3)
% set(h,'MarkerEdgeColor',clr)
% set(h,'MarkerFaceColor',clr)






function modPlot
set(gca,'TickDir','out')
set(gca,'TickLength',[.02 .02])
%set(gca,'YLim',[0 max(heights)])
%set(gca,'YTick',[0:1e-5:max(heights)])
%tmp=get(gca,'XLim');
%set(gca,'XLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
%tmp=get(gca,'YLim');
%set(gca,'YLim',[tmp(1)-range(tmp)/20 tmp(2)]); clear tmp
%axis square
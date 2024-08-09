%----------------Time-localized wavelet phase coherence--------------------
% TPC = tlphcoh(WT1,WT2,freq,fs,Optional:numcycles)
% calculates time-localized wavelet phase coherence TPC.
%
% Input:
% WT1,WT2 - wavelet transforms of two signals
% freq - frequencies used in wavelet transform
% fs - sampling frequency of a signals from which WT1, WT2 were calculated
% numcycles - number of cycles for calculating TPC (determines adaptive
%             window length, i.e. at 0.1 Hz it will be (1/0.1)*numcycles
%             seconds); default=10.
%
% Author: Dmytro Iatsenko (http://www.physics.lancs.ac.uk/research/nbmphysics/diats)
%--------------------------------------------------------------------------


function [TPC, phcoh] = tlphcoh_nina(TFR1,TFR2,freq,n,m,fs,plt, event_tbl, phcoh_surrogates, fig_title, visible, varargin)
%if plt=1 plot the TPC
[NF,L]=size(TFR1);
if nargin>12, wsize=varargin{1}; else wsize=10; end

IPC=exp(1i*angle(n*TFR1.*conj(m*TFR2)));
ZPC=IPC; ZPC(isnan(ZPC))=0; cumPC=[zeros(NF,1),cumsum(ZPC,2)];
TPC=zeros(NF,L)*NaN;
for fn=1:NF
    cs=IPC(fn,:); cumcs=cumPC(fn,:);
    tn1=find(~isnan(cs),1,'first'); tn2=find(~isnan(cs),1,'last');
    
    window=round((wsize/freq(fn))*fs); window=window+1-mod(window,2); hw=floor(window/2);
    % disp(window)
    if ~isempty(tn1+tn2) && window<=tn2-tn1
    locpc=abs(cumcs(tn1+window:tn2+1)-cumcs(tn1:tn2-window+1))/window;
    TPC(fn,tn1+hw:tn2-hw)=locpc;
    end
end
[phcoh,varargout] = wphcoh(n*TFR1,m*TFR2);


%%% plotting%%%%%%%

if plt==1

    close all;

    YY=freq; XX=(0:(L-1))/fs; ZZ=abs(TPC);  ZZ=ZZ.^2;

    scrsz=get(0,'ScreenSize'); figure('Position',[scrsz(3)/4,scrsz(4)/8,scrsz(3)/1.5,6*scrsz(4)/8], 'Visible', visible);

    axes('Position',[0.1,0.15,0.6,0.75],'Layer','top','YScale','log','Box','on','FontSize',16);
    TL=length(XX); FL=length(YY);
    pc=pcolor(XX,YY,ZZ); set(pc,'EdgeColor','none'); %title(ZZname);
    set(gca,'yscale','log')

    ylabel('Frequency (Hz)'); xlabel('Time (s)'); sgtitle(fig_title);
    xlim([0,(L-1)/fs]); ylim([freq(1),freq(end)+1]);
    hold on;
    ax = gca;
    % clim([0 ax.CLim(2)])
    %mycolormap = customcolormap(linspace(0,1,6), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd'});
    % mycolormap = customcolormap(linspace(0,1,6), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3'});
    mycolormap = jet;
    cb = colorbar('horiz');
    colormap(mycolormap);
    cb.Position = [cb.Position(1), cb.Position(2)-0.1, cb.Position(3), cb.Position(4)];

    axis tight
    ax = axis;
    sirnina_crt=0.5;
    frekvence = [7 15 30 50];
    y_tick_labels = {'1', '4', '7', '15', '30', '50', '100', '150'};
    y_ticks = [1 4 7 15 30 50 100 150];

    for ifreq = 1:size(frekvence, 2)
        l = yline(frekvence(ifreq));
        set(l,'color','black','linewidth',sirnina_crt,'linestyle','--');
    end
    yticklabels(y_tick_labels)
    yticks(y_ticks)

    % added by nina
    if ~isempty(event_tbl)
        for i = 1:size(event_tbl, 1)
            latency = cell2mat(event_tbl.latency(i)) / fs;
            l = line([latency latency], gca().YLim);
            set(l,'color','k','linewidth',sirnina_crt,'linestyle','--');
            if i == 1
                mp = latency;
            elseif rem(i, 2) == 0
                mp = latency-2000/fs;
            else
                mp = latency+2000/fs;
            end
            t = text(mp, max(ylim), event_tbl.type_new{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
            set(t,'color', 'k')
        end
    end

    %         axes('Position',[0.1,0.7,0.6,0.15],'Layer','top','Box','on','FontSize',16);
    %         plot(sgs1);axis tight;
    %         set(gca,'XTickLabel',{}); ylabel('Signal (t)');
    %         axes('Position',[0.1,0.7,0.6,0.25],'Layer','top','Box','on','FontSize',16);
    %         plot(sgs2);axis tight;
    %         set(gca,'XTickLabel',{}); ylabel('Signal (t)');

    %
    axes('Position',[0.75,0.15,0.2,0.75],'Layer','top','YScale','log','Box','on','FontSize',16);
    % mline=plot(phcoh',YY','-k','LineWidth',2); xlabel({'Time-averaged WPC'});
    set(gca,'yscale','log','YTickLabel', y_tick_labels, 'YTick', y_ticks)
    ylim([freq(1),freq(end)+1]);
    xlim([0,0.2]);

    %axis tight
    ax = axis;
    for ifreq = 1:size(frekvence, 2)
        l = yline(frekvence(ifreq));
        set(l,'color','black','linewidth',sirnina_crt,'linestyle','--');
    end
    yticklabels(y_tick_labels)
    yticks(y_ticks)

    if ~isempty(phcoh_surrogates)
        true_signal = phcoh';
        phcoh_surr_mean = mean(phcoh_surrogates, 1);
        phcoh_surr_2std = phcoh_surr_mean + 2*std(phcoh_surrogates, 1);
        n_phcohs = size(phcoh_surrogates, 1);
        hold on;
        % for idx_phcoh = 1: n_phcohs
        %     plot(phcoh_surrogates(idx_phcoh, :), YY', 'color', '#D3D3D3', 'LineStyle', '-');
        % end
        plot(phcoh_surr_mean, YY', 'color', 'b', 'LineStyle', '-', 'LineWidth', 2);
        % plot(phcoh_surr_2std, YY', 'color', '#A9A9A9', 'LineStyle', '--', 'LineWidth', 2)
        plot(true_signal, YY', 'k', 'LineWidth', 2)
        hold on;
        shade_yx(gca, YY', true_signal, YY', phcoh_surr_2std,'FillType',[1 2], 'FillColor', {'red'}, 'FillAlpha', 0.5)
    end
end

end


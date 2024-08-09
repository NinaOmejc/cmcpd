% plot single cwt results
fs=50;
channel_newname = 'avgcen';
subject='PDH09';
subject_group = subject(3);
task_split = 'SL';

load(['D:\Experiments\corticomuscular_analysis\data\real\cwt\IZO_v1_avgbrain_eqsplit' task_split '\' subject '\' subject '_IZO_eemg_cwt_' channel_newname '.mat'])
cwt_avg = WT;
freqs = f;

fig_title = ['cwts avg | ' task_split ' | ' char(subject) ' | ' channel_newname];
freq_mask = (freqs>=4) & (freqs <= 100);
freq = freqs(freq_mask);
YY = freq;
L = length(cwt_avg);
XX = (0:(L-1))/fs; 
ZZ = cwt_avg(freq_mask, :);

scrsz=get(0,'ScreenSize'); 
figure('Position',[scrsz(3)/4,scrsz(4)/8,scrsz(3)/1.5,6*scrsz(4)/8], 'Visible','on');

axes('Position',[0.1,0.15,0.6,0.53],'Layer','top','YScale','log','Box','on','FontSize',16);
TL=length(XX); FL=length(YY);
pc=pcolor(XX,YY,ZZ); 
set(pc,'EdgeColor','none'); %title(ZZname);
ylabel('Frequency (Hz)'); xlabel('Time (s)');
xlim([0,(L-1)/fs]); ylim([freq(1),freq(end)]);
set(gca,'yscale','log')

ax = gca;
% mycolormap = customcolormap(linspace(0,1,6), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd'});
% mycolormap = customcolormap(linspace(0,1,6), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3'});
mycolormap = jet;
cb = colorbar('horiz');
colormap(mycolormap);
cb.Position = [cb.Position(1), cb.Position(2)-0.1, cb.Position(3), cb.Position(4)];
clim([0 0.01*ax.CLim(2)])

hold on; axis tight
ax = axis;
sirnina_crt=0.5;
frekvence = [4 7.5 15 30 50 60 100];
for ifreq = 1:size(frekvence, 2)
    l = yline(frekvence(ifreq));
    set(l,'color','black','linewidth',sirnina_crt,'linestyle','--');
end
yticklabels({'4', '7.5', '15', '30', '50', '60', '100'})
yticks(frekvence)

% plot average
TL=length(XX); FL=length(YY);
mx=zeros(FL,1); for fn=1:FL, mx(fn)=mean(ZZ(fn,~isnan(ZZ(fn,:))),2); end
ZZname = 'WT';
axes('Position',[0.75,0.15,0.2,0.53],'Layer','top','YScale','log','Box','on','FontSize',16);
mline=plot(mx',YY','-k','LineWidth',2); xlabel({'Time-averaged', ZZname});
set(gca,'yscale','log','YTickLabel',{})
ylim([freq(1),freq(end)]);
sgtitle(fig_title);

axis tight
ax = axis;
for ifreq = 1:size(frekvence, 2)
    l = line([ax(1) ax(2)],[frekvence(ifreq) frekvence(ifreq)]);
    set(l,'color','k','linewidth',sirnina_crt,'linestyle','--');
end
yticklabels({'4', '7.5', '15', '30', '50', '60', '100'})
yticks(frekvence)


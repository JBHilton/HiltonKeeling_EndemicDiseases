% This script makes heatmaps of contact behaviour. The commented out
% section is for the "all contacts" matrices.

load('Parameters/UK_MixingData.mat')

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
imagesc(24*D_All);
set(gca,'YDir','normal');
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
yticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
% c=colorbar;
% c.Label.String='Daily exposure (hours)';
% c.AxisLocation='in';
% c.Position(1)=c.Position(1)+0.07;
% caxis([0 15]);
% c.Ticks=0:3:15;
xlabel('Participant age class');
ylabel('Contact age class');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
imagesc(24*D_Ext');
set(gca,'YDir','normal');
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
yticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
caxis([0 18]);
% c.Ticks=0:3:15;
xlabel('Participant age class');
ylabel('Contact age class');

load('Parameters/Kenya_MixingData.mat');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
imagesc(24*D_All);
set(gca,'YDir','normal');
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
yticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
caxis([0 18]);
xlabel('Participant age class');
ylabel('Contact age class');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
i=imagesc(24*D_Ext');
set(gca,'YDir','normal');
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
yticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Participant age class');
ylabel('Contact age class');
c=colorbar;
c.Label.String='Daily exposure (hours)';
c.AxisLocation='in';
caxis([0 18]);
c.Ticks=0:6:18;
ax=gca;
c.Position(1)=c.Position(1)+0.07;
c.Position(2)=ax.Position(2);
c.Position(4)=ax.Position(4)-0.02;
% This script plots the results of our sensitivity analysis on the
% structure parameter sigma.

range=1:11;

%% First do mumps in the UK

load('ModelOutput/UKMumpsSigmaLoop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Prev_AgeOnly(range),'-sg',sigma(range),Prev_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 1.7e-4 2e-4])
axis square
xlabel('Structure parameter \sigma')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
prevlabel=text(0.5,0.5,'$\ \bar{I}$','Interpreter','latex','FontSize',24);
prevlabel.Rotation=90;
prevlabel.VerticalAlignment='baseline';
yl.VerticalAlignment='baseline';
yl.HorizontalAlignment='right';
prevlabel.HorizontalAlignment='left';
pause(0.5);
prevlabel.Position(2)=yl.Position(2);
pause(0.5);
prevlabel.Position(1)=yl.Position(1);
vdiff=prevlabel.Extent(2)-yl.Extent(2);
hdiff=prevlabel.Extent(1)-yl.Extent(1);
delete(yl);
yl=ylabel('Equilibrium prevalence');
prevlabel.HorizontalAlignment='left';
prevlabel.VerticalAlignment='top';
prevlabel.Position(1)=yl.Extent(1)+hdiff;
prevlabel.Position(2)=yl.Extent(2)+vdiff;

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),r_AgeOnly(range),'-sg',sigma(range),r_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 0 0.5])
axis square
xlabel('Structure parameter \sigma')
ylabel('Growth rate r')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Rs_AgeOnly(range),'-sg',sigma(range),Rs_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 2 3])
axis square
xlabel('Structure parameter \sigma')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
l=legend('No Household Structure','Household Structure','Location','southeast');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),PR_AgeOnly(range),'-sg',sigma(range),PR_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 0 1])
axis square
xlabel('Structure parameter \sigma')
ylabel('Childhood infection probability P_R')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

%% Now do measles in the UK

load('ModelOutput/UKMeaslesSigmaLoop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Prev_AgeOnly(range),'-sg',sigma(range),Prev_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 1.7e-4 2e-4])
axis square
xlabel('Structure parameter \sigma')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
prevlabel=text(0.5,0.5,'$\ \bar{I}$','Interpreter','latex','FontSize',24);
prevlabel.Rotation=90;
prevlabel.VerticalAlignment='baseline';
yl.VerticalAlignment='baseline';
yl.HorizontalAlignment='right';
prevlabel.HorizontalAlignment='left';
pause(0.5);
prevlabel.Position(2)=yl.Position(2);
pause(0.5);
prevlabel.Position(1)=yl.Position(1);
vdiff=prevlabel.Extent(2)-yl.Extent(2);
hdiff=prevlabel.Extent(1)-yl.Extent(1);
delete(yl);
yl=ylabel('Equilibrium prevalence');
prevlabel.HorizontalAlignment='left';
prevlabel.VerticalAlignment='top';
prevlabel.Position(1)=yl.Extent(1)+hdiff;
prevlabel.Position(2)=yl.Extent(2)+vdiff;

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),r_AgeOnly(range),'-sg',sigma(range),r_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 0 2])
axis square
xlabel('Structure parameter \sigma')
ylabel('Growth rate r')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Rs_AgeOnly(range),'-sg',sigma(range),Rs_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 9 16])
axis square
xlabel('Structure parameter \sigma')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),PR_AgeOnly(range),'-sg',sigma(range),PR_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 0 1])
axis square
xlabel('Structure parameter \sigma')
ylabel('Childhood infection probability P_R')
titletext=text(0.5,0.75,'UK, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
figname='PR_UKSigmaMeasles';

%% Now do mumps in Kenya

load('ModelOutput/KenyaMumpsSigmaLoop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Prev_AgeOnly(range),'-sg',sigma(range),Prev_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 2.9e-4 3.2e-4])
axis square
xlabel('Structure parameter \sigma')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
prevlabel=text(0.5,0.5,'$\ \bar{I}$','Interpreter','latex','FontSize',24);
prevlabel.Rotation=90;
prevlabel.VerticalAlignment='baseline';
yl.VerticalAlignment='baseline';
yl.HorizontalAlignment='right';
prevlabel.HorizontalAlignment='left';
pause(0.5);
prevlabel.Position(2)=yl.Position(2);
pause(0.5);
prevlabel.Position(1)=yl.Position(1);
vdiff=prevlabel.Extent(2)-yl.Extent(2);
hdiff=prevlabel.Extent(1)-yl.Extent(1);
delete(yl);
yl=ylabel('Equilibrium prevalence');
prevlabel.HorizontalAlignment='left';
prevlabel.VerticalAlignment='top';
prevlabel.Position(1)=yl.Extent(1)+hdiff;
prevlabel.Position(2)=yl.Extent(2)+vdiff;

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),r_AgeOnly(range),'-sg',sigma(range),r_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 0 0.5])
axis square
xlabel('Structure parameter \sigma')
ylabel('Growth rate r')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Rs_AgeOnly(range),'-sg',sigma(range),Rs_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 3 8])
axis square
xlabel('Structure parameter \sigma')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),PR_AgeOnly(range),'-sg',sigma(range),PR_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 0 1])
axis square
xlabel('Structure parameter \sigma')
ylabel('Childhood infection probability P_R')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

%% Now do measles in Kenya

load('ModelOutput/KenyaMeaslesSigmaLoop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Prev_AgeOnly(range),'-sg',sigma(range),Prev_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 2.9e-4 3.2e-4])
axis square
xlabel('Structure parameter \sigma')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
prevlabel=text(0.5,0.5,'$\ \bar{I}$','Interpreter','latex','FontSize',24);
prevlabel.Rotation=90;
prevlabel.VerticalAlignment='baseline';
yl.VerticalAlignment='baseline';
yl.HorizontalAlignment='right';
prevlabel.HorizontalAlignment='left';
pause(0.5);
prevlabel.Position(2)=yl.Position(2);
pause(0.5);
prevlabel.Position(1)=yl.Position(1);
vdiff=prevlabel.Extent(2)-yl.Extent(2);
hdiff=prevlabel.Extent(1)-yl.Extent(1);
delete(yl);
yl=ylabel('Equilibrium prevalence');
prevlabel.HorizontalAlignment='left';
prevlabel.VerticalAlignment='top';
prevlabel.Position(1)=yl.Extent(1)+hdiff;
prevlabel.Position(2)=yl.Extent(2)+vdiff;
l=legend('No Household Structure','Household Structure','Location','southeast');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),r_AgeOnly(range),'-sg',sigma(range),r_Full(range),'-x','LineWidth',1.5,'MarkerSize',10);
axis([0 1 0 2])
axis square
xlabel('Structure parameter \sigma')
ylabel('Growth rate r')
titletext=text(0.5,1.5,'Kenya, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
l=legend('No Household Structure','Household Structure','Location','southeast');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),Rs_AgeOnly(range),'-sg',sigma(range),Rs_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 10 50])
axis square
xlabel('Structure parameter \sigma')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim'); titletext=text(0.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(sigma(range),PR_AgeOnly(range),'-sg',sigma(range),PR_Full(range),'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 1 0 1])
axis square
xlabel('Structure parameter \sigma')
ylabel('Childhood infection probability P_R')
titletext=text(0.5,0.8,'Kenya, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
l=legend('No Household Structure','Household Structure','Location','southeast');
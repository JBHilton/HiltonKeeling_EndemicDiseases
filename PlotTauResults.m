% This script plots the results of our sensitivity analysis on the
% transmission rate tau.

%% First do gamma=7 in the UK

load('ModelOutput/UKTau7Loop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Prev_Rand,'-or',tau,Prev_POLY,'-sg',tau,Prev_HH,'-^m',tau,Prev_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 4.5e-4])
axis square
xlabel('Transmission rate \tau')
yl=ylabel('Equilibrium prevalence');
pause(1)
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=7','FontSize',24);
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
plot(tau,r_Rand,'-or',tau,r_POLY,'-sg',tau,r_HH,'-^m',tau,r_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 -0.2 2.5])
axis square
pause(1)
xlabel('Transmission rate \tau')
ylabel('Growth rate r')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=7','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Rs_Rand,'-or',tau,Rs_POLY,'-sg',tau,Rs_HH,'-^m',tau,Rs_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis square
pause(1)
xlabel('Transmission rate \tau')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=7','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,PR_Rand,'-or',tau,PR_POLY,'-sg',tau,PR_HH,'-^m',tau,PR_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 1])
pause(1)
axis square
xlabel('Transmission rate \tau')
ylabel('Childhood infection probability P_R')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=7','FontSize',24);
titletext.Position(1)=titletext.Position(1)-1.1*titletext.Extent(3);
titletext.Position(2)=titletext.Position(2)-0.5*titletext.Extent(4);

%% Now do gamma=14 in the UK
clear
load('ModelOutput/UKTau14Loop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Prev_Rand,'-or',tau,Prev_POLY,'-sg',tau,Prev_HH,'-^m',tau,Prev_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 4.5e-4])
axis square
xlabel('Transmission rate \tau')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=14','FontSize',24);
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
plot(tau,r_Rand,'-or',tau,r_POLY,'-sg',tau,r_HH,'-^m',tau,r_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 -0.2 2.5])
axis square
xlabel('Transmission rate \tau')
ylabel('Growth rate r')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=14','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Rs_Rand,'-or',tau,Rs_POLY,'-sg',tau,Rs_HH,'-^m',tau,Rs_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis square
xlabel('Transmission rate \tau')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'UK, \gamma^{-1}=14','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,PR_Rand,'-or',tau,PR_POLY,'-sg',tau,PR_HH,'-^m',tau,PR_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 1])
axis square
pause(1)
xlabel('Transmission rate \tau')
ylabel('Childhood infection probability P_R')
titletext=text(1,0.8,'UK, \gamma^{-1}=14','FontSize',24);

%% Now do gamma=7 in Kenya
clear
load('ModelOutput/KenyaTau7Loop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Prev_Rand,'-or',tau,Prev_POLY,'-sg',tau,Prev_HH,'-^m',tau,Prev_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 7e-4])
axis square
xlabel('Transmission rate \tau')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'Kenya, \gamma^{-1}=7','FontSize',24);
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
plot(tau,r_Rand,'-or',tau,r_POLY,'-sg',tau,r_HH,'-^m',tau,r_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 -0.2 2.5])
axis square
xlabel('Transmission rate \tau')
ylabel('Growth rate r')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'Kenya, \gamma^{-1}=7','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Rs_Rand,'-or',tau,Rs_POLY,'-sg',tau,Rs_HH,'-^m',tau,Rs_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis square
pause(1)
xlabel('Transmission rate \tau')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'Kenya, \gamma^{-1}=7','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
l=legend('Homogeneous Mixing','Age-Structured','Household-Structured','Full Model','Location','west');
l.Position(2)=1.1*l.Position(2);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,PR_Rand,'-or',tau,PR_POLY,'-sg',tau,PR_HH,'-^m',tau,PR_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 1])
axis square
pause(1)
xlabel('Transmission rate \tau')
ylabel('Childhood infection probability P_R')
titletext=text(1,0.8,'Kenya, \gamma^{-1}=7','FontSize',24);

%% Now do gamma=14 in Kenya
clear
load('ModelOutput/KenyaTau14Loop.mat');
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Prev_Rand,'-or',tau,Prev_POLY,'-sg',tau,Prev_HH,'-^m',tau,Prev_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 7e-4])
axis square
xlabel('Transmission rate \tau')
yl=ylabel('Equilibrium prevalence');
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'Kenya, \gamma^{-1}=14','FontSize',24);
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
l=legend('Homogeneous Mixing','Age-Structured','Household-Structured','Full Model','Location','southeast');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,r_Rand,'-or',tau,r_POLY,'-sg',tau,r_HH,'-^m',tau,r_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 -0.2 2.5])
axis square
xlabel('Transmission rate \tau')
ylabel('Growth rate r')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'Kenya, \gamma^{-1}=14','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
l=legend('Homogeneous Mixing','Age-Structured','Household-Structured','Full Model','Location','southeast');

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,Rs_Rand,'-or',tau,Rs_POLY,'-sg',tau,Rs_HH,'-^m',tau,Rs_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis square
xlabel('Transmission rate \tau')
ylabel('Household reproduction ratio R^*')
ylim=get(gca,'ylim');
titletext=text(1,0.99*ylim(2),'Kenya, \gamma^{-1}=14','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
plot(tau,PR_Rand,'-or',tau,PR_POLY,'-sg',tau,PR_HH,'-^m',tau,PR_Full,'-x','LineWidth',1.5,'MarkerSize',10)
axis([0 2 0 1])
axis square
pause(1)
xlabel('Transmission rate \tau')
ylabel('Childhood infection probability P_R')
ylim=get(gca,'ylim');
titletext=text(1,0.8,'Kenya, \gamma^{-1}=14','FontSize',24);
l=legend('Homogeneous Mixing','Age-Structured','Household-Structured','Full Model','Location','southeast');
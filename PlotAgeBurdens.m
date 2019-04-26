% This script generates histograms of infectious prevalence by age
% class, for the UK under the four transmission structures and for the four
% disease-demography combinations.

%% Start with UK histograms

load('Parameters/UK_MixingData.mat');

load('Parameters/UKDemography.mat');

load('ModelOutput/UKStructureEqDists.mat');

load('ModelOutput/UKMeaslesEqDist.mat');

load('ModelOutput/UKMumpsEqDist.mat');

nVectI=nVect(2,:);
% The cells in Index are the vectors locating the demographic states
Index=cell(1,length(NGrid));
I_T_mumps_UK=zeros(1,length(NGrid));
I_T_measles_UK=zeros(1,length(NGrid));
I_T_Hom=zeros(1,length(NGrid));
I_T_POLY=zeros(1,length(NGrid));
I_T_HH=zeros(1,length(NGrid));
I_T_Full=zeros(1,length(NGrid));
for i=1:length(NGrid)
    Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
    I_T_mumps_UK(i)=nVectI(Index{i})*Equil_MumpsUK(Index{i})/(H_T(i)*NGrid(i));
    I_T_measles_UK(i)=nVectI(Index{i})*Equil_MeaslesUK(Index{i})/(H_T(i)*NGrid(i));
    I_T_Hom(i)=nVectI(Index{i})*Equil_Hom(Index{i})/(H_T(i)*NGrid(i));
    I_T_POLY(i)=nVectI(Index{i})*Equil_POLY(Index{i})/(H_T(i)*NGrid(i));
    I_T_HH(i)=nVectI(Index{i})*Equil_HH(Index{i})/(H_T(i)*NGrid(i));
    I_T_Full(i)=nVectI(Index{i})*Equil_Full(Index{i})/(H_T(i)*NGrid(i));
end

PE_Sum = H_T*E;
IC_mumps=100*((H_T.*I_T_mumps_UK)*E)./PE_Sum;
IC_measles=100*((H_T.*I_T_measles_UK)*E)./PE_Sum;
IC_Hom=100*((H_T.*I_T_Hom)*E)./PE_Sum;
IC_POLY=100*((H_T.*I_T_POLY)*E)./PE_Sum;
IC_HH=100*((H_T.*I_T_HH)*E)./PE_Sum;
IC_Full=100*((H_T.*I_T_Full)*E)./PE_Sum;

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_mumps);
axis([0.4 6.6 0 0.12]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_measles);
axis([0.4 6.6 0 0.12]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'UK, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_Hom);
axis([0.4 6.6 0 0.04]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Homogeneous Mixing','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_POLY);
axis([0.4 6.6 0 0.04]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Age Structured Mixing','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_HH);
axis([0.4 6.6 0 0.04]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Household-Structured Mixing','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_Full);
axis([0.4 6.6 0 0.04]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Full Model','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

%% Now do Kenya histograms

load('Parameters/Kenya_MixingData.mat');

load('Parameters/KenyaDemography.mat');

load('ModelOutput/KenyaMeaslesEqDist.mat');

load('ModelOutput/KenyaMumpsEqDist.mat');

nVectI=nVect(2,:);
% The cells in Index are the vectors locating the demographic states
Index=cell(1,length(NGrid));
I_T_mumps_Kenya=zeros(1,length(NGrid));
I_T_measles_Kenya=zeros(1,length(NGrid));
for i=1:length(NGrid)
    Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
    I_T_mumps_Kenya(i)=nVectI(Index{i})*Equil_MumpsKenya(Index{i})/(H_T(i)*NGrid(i));
    I_T_measles_Kenya(i)=nVectI(Index{i})*Equil_MeaslesKenya(Index{i})/(H_T(i)*NGrid(i));
end

PE_Sum = H_T*E;
IC_mumps=100*((H_T.*I_T_mumps_Kenya)*E)./PE_Sum;
IC_measles=100*((H_T.*I_T_measles_Kenya)*E)./PE_Sum;

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_mumps);
axis([0.4 6.6 0 0.12]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Mumps','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);

fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
bar(1:6,IC_measles);
axis([0.4 6.6 0 0.12]);
axis square;
xticklabels({'C_1', 'C_2', 'C_3', 'C_4', 'C_5', 'C_6'});
xlabel('Age class');
ylabel('Prevalence (%)');
ylim=get(gca,'ylim'); titletext=text(3.5,ylim(1)+0.99*(ylim(2)-ylim(1)),'Kenya, Measles','FontSize',24);
titletext.Position(1:2)=titletext.Position(1:2)-0.5*titletext.Extent(3:4);
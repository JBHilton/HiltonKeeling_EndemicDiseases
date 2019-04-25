% Calculates disease-free equilibria for UK and Kenya and plots size/age
% distribution as stacked bar charts.

% Define the set of colours we will use to draw stratified prevalence plots
ColourMap=jet;


% First do UK
load('Parameters\UK_MixingData.mat') % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

% Define demographic parameters
TB=3*365;
kB=2;
TL=25*365;
kL=4;
TD=80*365;
kR=5;
Distrib_Children=[0.17 0.18 0.37 0.17 0.1];
Distrib_Children=Distrib_Children/sum(Distrib_Children);
Exp_Children = (0:length(Distrib_Children)-1)*Distrib_Children'; % Expected number of children
TR = TD-TL-TL-TB*Exp_Children; % Reset rate = age of death - expected time with kids
tmp=[0 cumsum(Distrib_Children)];
StopProb=[0 Distrib_Children./(1-tmp(1:(end-1)))]; % Prob of stopping after each number of kids
StopProb(end)=1;
maxN = find(StopProb==1,1); % Can't have more kids after this
if isempty(maxN)
    StopProb(end)=1;
    maxN=length(StopProb);
end

% We now calculate the disease-free equilibrium to get distribution of
% household sizes
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(0,0, 0, 0, 1, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0);
        flag = 1;
    catch
        warning('error finding disease-free equilibrium')
    end
end
H_DiseaseFree=H_Eq;
nVectN=sum(nVect,1);
% The cells in Index are the vectors locating the states we are stratifying
% by
Index=cell(1,2*maxN-3);
x_lab_UK=cell(1,2*maxN-3);
Index{1}=find(nVectN==2&nTicker<kB);
x_lab_UK{1}=num2str(2);
Index{2*maxN-3}=find(nVectN==2&nTicker>=kB);
x_lab_UK{2*maxN-3}=num2str(2);
Index{maxN-1}=find(nVectN==maxN);
x_lab_UK{maxN-1}=num2str(maxN);
for i=3:maxN-1
    Index{i-1}=find(nVectN==i&nTicker<kB+kL);
    x_lab_UK{i-1}=num2str(i);
    Index{2*(maxN)-i-1}=find(nVectN==i&nTicker>=kB+kL);
    x_lab_UK{2*(maxN)-i-1}=num2str(i);
end
y_UK=zeros(1,2*maxN-3);
for i=1:2*maxN-3
    y_UK(i)=100*sum(H_DiseaseFree(Index{i}));
end

% Now do the plot
UKHHDistList=figure('DefaultAxesFontSize',24);
b_young=bar(1:maxN-1,y_UK(1:maxN-1));
hold on
b_old=bar(maxN:2*maxN-3,y_UK(maxN:end));
set(b_young,'FaceColor',ColourMap(1,:));
set(b_old,'FaceColor',ColourMap(end,:));
hold on
legend('k<k_B+k_L','k\geq k_B+k_L','Location','northeast')
line([maxN-1,maxN-1],[0,45],'Color','k','LineStyle','--','LineWidth',1) % Divider between young and old households
axis([0 2*maxN-2 0 45])
axis square
xticks(1:2*maxN-3)
xticklabels(x_lab_UK)
xlabel('Household size')
ylabel('Frequency (%)')


%%
% Now do Kenya
load('Parameters\Kenya_MixingData.mat') % Contains d_ext, d_int, ClassProb, D_Alld, D_Ext, E, NGrid, tickGrid, DemGrid

% Define demographic parameters
TB=(35.2/12)*365;
kB=2;
TL=22.7*365;
kL=4;
TD=67*365;
kR=5;
Distrib_Children=[0.019 0.038 0.098 0.129 0.15 0.132 0.119 0.108 0.078 0.05 0.079];
Constrained_max=floor(TL/TB)+1; % This is the most kids we can have under our model assumptions with these parameter choices
Distrib_Children(Constrained_max+1)=sum(Distrib_Children(Constrained_max+1:end)); % +1 since first one is no kids
Distrib_Children=Distrib_Children(1:Constrained_max+1);
Distrib_Children=Distrib_Children/sum(Distrib_Children);
Exp_Children = (0:length(Distrib_Children)-1)*Distrib_Children'; % Expected number of children
TR = TD-TL-TL-TB*Exp_Children; % Reset rate = age of death - expected time with kids
tmp=[0 cumsum(Distrib_Children)];
StopProb=[0 Distrib_Children./(1-tmp(1:(end-1)))]; % Prob of stopping after each number of kids
StopProb(end)=1;
maxN = find(StopProb==1,1); % Can't have more kids after this
if isempty(maxN)
    StopProb(end)=1;
    maxN=length(StopProb);
end

% We now calculate the disease-free equilibrium to get distribution of
% household sizes
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(0,0, 0, 0, 1, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0);
        flag = 1;
    catch
        warning('error finding disease-free equilibrium')
    end
end
H_DiseaseFree=H_Eq;
nVectN=sum(nVect,1);
nVectI=nVect(2,:);
% The cells in Index are the vectors locating the states we are stratifying
% by
Index=cell(1,2*maxN-3);
x_lab_Kenya=cell(1,2*maxN-3);
Index{1}=find(nVectN==2&nTicker<kB);
x_lab_Kenya{1}=num2str(2);
Index{2*maxN-3}=find(nVectN==2&nTicker>=kB);
x_lab_Kenya{2*maxN-3}=num2str(2);
Index{maxN-1}=find(nVectN==maxN);
x_lab_Kenya{maxN-1}=num2str(maxN);
for i=3:maxN-1
    Index{i-1}=find(nVectN==i&nTicker<kB+kL);
    x_lab_Kenya{i-1}=num2str(i);
    Index{2*(maxN)-i-1}=find(nVectN==i&nTicker>=kB+kL);
    x_lab_Kenya{2*(maxN)-i-1}=num2str(i);
end
y_Kenya=zeros(1,2*maxN-3);
for i=1:2*maxN-3
    y_Kenya(i)=100*sum(H_DiseaseFree(Index{i}));
end

KenyaHHDistList=figure('DefaultAxesFontSize',24);
b_young=bar(1:maxN-1,y_Kenya(1:maxN-1));
hold on
b_old=bar(maxN:2*maxN-3,y_Kenya(maxN:end));
set(b_young,'FaceColor',ColourMap(1,:));
set(b_old,'FaceColor',ColourMap(end,:));
hold on
legend('k<k_B+k_L','k\geq k_B+k_L','Location','northeast')
line([maxN-1,maxN-1],[0,45],'Color','k','LineStyle','--','LineWidth',1) % Divider between young and old households
axis([0 2*maxN-2 0 45])
axis square
xticks(1:2*maxN-3)
xticklabels(x_lab_Kenya)
xlabel('Household size')
ylabel('Frequency (%)')
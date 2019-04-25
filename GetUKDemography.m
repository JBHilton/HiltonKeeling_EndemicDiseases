% In this script we define all of the demographic parameters for the UK and
% calculate the equilibrium distribution of demographic states and the
% matrices E and \tilde{E} (coded as E1 and E2 respectively). Having this
% data saved means we don't have to recalculate it at the beginning of all
% the other scripts.

load('Parameters/UK_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

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
        [~,~,H_T,~,H_Eq, nVect, nTicker, ~, ~, ~]=HH_demo_structured(0,0, 0, 0, 1, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0);
        flag = 1;
    catch me
        warning('error finding disease-free equilibrium')
    end
end
DiseaseFree=H_Eq;
nVectN=sum(nVect,1);

% The following loop makes sure that all the rows of E sum to the household
% size. This is necessary because we calculated it using Monte Carlo
% integration, which will not be completely accurate.
for i=1:length(NGrid)
    if NGrid(i)>2
        E(i,end)=E(i,end)-2;
        E(i,:)=(NGrid(i)-2)*E(i,:)/sum(E(i,:));
        E(i,end)=E(i,end)+2;
    end
end
PQ_Sum = H_T*E;
E2=H_T.*E'./PQ_Sum';

filename=['Parameters/UKDemography_' datestr(now,'ddmmyy_HHMMSS')];
save(filename, 'TB', 'kB', 'TL', 'kL', 'TD', 'kR', 'Exp_Children', 'TR',...
    'StopProb', 'maxN', 'H_T', 'DiseaseFree', 'nVect','nTicker',...
    'nVectN','E','E2');
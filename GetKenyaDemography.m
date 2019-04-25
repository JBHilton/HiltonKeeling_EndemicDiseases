% In this script we define all of the demographic parameters for the UK and
% calculate the equilibrium distribution of demographic states and the
% matrices E and \tilde{E} (coded as E1 and E2 respectively). Having this
% data saved means we don't have to recalculate it at the beginning of all
% the other scripts.

load('Parameters/Kenya_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

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
        [~,~,H_T,~,H_Eq, nVect, nTicker, ~, ~, ~]=HH_demo_structured(0,0, 0, 0, 1, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0);
        flag = 1;
    catch meQ
        warning('error finding disease-free equilibrium')
        %rethrow(me);
    end
end
DiseaseFree=H_Eq;
nVectN=sum(nVect,1);

% The following loop makes sure that all the rows of E sum to the household
% size. This is necessary because we calculated it using Monte Carlo
% integration, which will not be completely accurate.
for i=1:length(NGrid)
    if NGrid(i)>2
        E(i,end-1)=E(i,end-1)-2; % Because of elder class, parents are in second-to-one class
        E(i,:)=(NGrid(i)-2)*E(i,:)/sum(E(i,:));
        E(i,end)=E(i,end)+2;
    end
end
PQ_Sum = H_T*E;
E2=H_T.*E'./PQ_Sum';

filename=['Parameters/KenyaDemography_' datestr(now,'ddmmyy_HHMMSS')];
save(filename, 'TB', 'kB', 'TL', 'kL', 'TD', 'kR', 'Exp_Children', 'TR',...
    'StopProb', 'maxN', 'H_T', 'DiseaseFree', 'nVect','nTicker',...
    'nVectN','E','E2');
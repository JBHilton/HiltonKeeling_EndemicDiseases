% Calculates endemic equilibria for MUMPS in UK with full age- and
% household-structured model, then calculates it for other three
% combinations of age- and household-structure with equilbalent R0.
% The outputs of this script are plotted in Figure 3 of the paper and
% Figure 1 of the supplement, and listed in Table 2 of the paper.

load('Parameters/UK_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/UKDemography.mat')

% Start with full age- and household-structured model
g=1/8; % Infectious period
d_int=sum(ClassProb*(D_All-D_Ext));
d_ext=sum(ClassProb*D_Ext);
tau=-log(0.689)*g/d_int_kids; % Unit time transmission rate based on Hope-Simpson
beta_int=tau*d_int; % Internal contact rate

% Do scaling of WAIFW matrix by multiplying through by transmission rate
Inf_Ext=tau*D_Ext; % Rates of external infection

% Start by calculating early growth parameters to fit to later
flag = 0;
while flag==0
    try
        [Rs_Full, r_Full] = Grab_Epi_Data(E,E2,Inf_Ext,0,1, 0, beta_int, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
        flag = 1;
    catch
        warning('Error in Grab_Epi_Data.m for unvaccinated population');
    end
end

R0_Full=1+r_Full/g;

% Now do equilibrium prevalence
iter = 1; % Helps us keep track of how far along we are, if needed
P_R=0;
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_int, Q_ext, Q_demo]=HH_demo_structured(1e-6,1e-6, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
        flag = 1;
    catch
        warning('error finding initial guess state')
    end
end

% Now do recursive process to find equilibrium
flag = 0;
while flag==0
    % Note that we have no unstructured infection here
    Struc_Inf = E*Inf_Ext*E2*I_T'; % Age-structured ext. inf.
    try
        [I_bar_new,I_T_new,~,P_R_new,H_Eq_new, ~, ~, ~, ~, ~]=HH_demo_structured(0,Struc_Inf, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
        if sum(abs(H_Eq-H_Eq_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_T_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-P_R_new)/P_R_new<1e-6
            flag=1;
        end
        I_bar = I_bar_new; I_T=I_T_new; P_R = P_R_new; H_Eq = H_Eq_new; % Get new ambient conditions
    catch
        warning(['error on iteration' num2str(iter)]);
    end
    iter = iter+1;
    if iter>100
        flag = 1;
        warning('Failed to converge for i=1');
    end
end
Prev_Full=I_bar; PR_Full=P_R;
Equil_Full=H_Eq; % Equilibrium distribution for full age- and household-structured model

%% Now do homogeneous mixing

beta_Hom=g*R0_Full; % Unstructured external contact rate
tau_Hom=beta_Hom/(d_int+d_ext);

flag = 0;
while flag==0
    try
        [Rs_Hom, r_Hom] = Grab_Epi_Data(E,E2,Inf_Ext,beta_Hom,0, 0, 0, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
        flag = 1;
    catch
        warning('Error in Grab_Epi_Data.m for unvaccinated population');
    end
end
R0_Hom=1+r_Hom/g; % Should be same as R0_Full

% Now do equilibrium prevalence
iter = 1; % Helps us keep track of how far along we are, if needed
P_R=0;
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_,P_R,H_Eq, nVect, nTicker, Q_int, Q_ext, Q_demo]=HH_demo_structured(1e-6,1e-6, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
        flag = 1;
    catch
        warning('error finding initial guess state')
    end
end

% Now do recursive process to find equilibrium
flag = 0;
while flag==0
    Unstr_Inf = I_bar*beta_Hom*ones(1,maxN-1);
    try
        [I_bar_new,I_T_new,~,P_R_new,H_Eq_new, ~, ~, ~, ~, ~]=HH_demo_structured(Unstr_Inf,0, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
        if sum(abs(H_Eq-H_Eq_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_T_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-P_R_new)/P_R_new<1e-6
            flag=1;
        end
        I_bar = I_bar_new; I_T=I_T_new; P_R = P_R_new; H_Eq = H_Eq_new; % Get new ambient conditions
    catch
        warning(['error on iteration' num2str(iter)]);
    end
    iter = iter+1;
    if iter>100
        flag = 1;
        warning('Failed to converge for i=1');
    end
end
Prev_Hom=I_bar; PR_Hom=P_R;
Equil_Hom=H_Eq;

%% Now do POLYMOD age-structure without household structure

% Do scaling of WAIFW matrix by multiplying through by transmission rate
[Rs_POLY_high, r_POLY_high] = Grab_Epi_Data(E,E2,tau*D_All,0,1, 0, 0, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
Inf_All=(tau-tau*(1+r_POLY_high/g-R0_Full)/(1+r_POLY_high/g))*D_All; % R0 is linear in tau, this finds the middle assuming R0(tau=0)=0
tau_POLY=(tau-tau*(1+r_POLY_high/g-R0_Full)/(1+r_POLY_high/g));

flag = 0;
while flag==0
    try
        [Rs_POLY, r_POLY] = Grab_Epi_Data(E,E2,Inf_All,0,1, 0, 0, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
        flag = 1;
    catch
        warning('Error in Grab_Epi_Data.m for unvaccinated population');
    end
end
R0_POLY=1+r_POLY/g; % Should be same as R0_Full

% Now do equilibrium prevalence
iter = 1; % Helps us keep track of how far along we are, if needed
P_R=0;
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_,P_R,H_Eq, nVect, nTicker, Q_int, Q_ext, Q_demo]=HH_demo_structured(1e-6,1e-6, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
        flag = 1;
    catch
        warning('error finding initial guess state')
    end
end

% Now do recursive process to find equilibrium
flag = 0;
while flag==0
        Struc_Inf = E*Inf_All*E2*I_T'; % Age-structured ext. inf.
    try
        [I_bar_new,I_T_new,~,P_R_new,H_Eq_new, ~, ~, ~, ~, ~]=HH_demo_structured(0,Struc_Inf, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
        if sum(abs(H_Eq-H_Eq_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_T_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-P_R_new)/P_R_new<1e-6
            flag=1;
        end
        I_bar = I_bar_new; I_T=I_T_new; P_R = P_R_new; H_Eq = H_Eq_new; % Get new ambient conditions
    catch
        warning(['error on iteration' num2str(iter)]);
    end
    iter = iter+1;
    if iter>100
        flag = 1;
        warning('Failed to converge for i=1');
    end
end
Prev_POLY=I_bar; PR_POLY=P_R;
Equil_POLY=H_Eq;

%% Now do household structure with no age structure

% Use interval bisection to find appropriate tau

beta_int=tau*d_int;
beta_ext=tau*d_ext;
[Rs_HH_low, r_HH_low] = Grab_Epi_Data(E,E2,0,beta_ext,0, 0, beta_int, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
R0_HH_low=1+r_HH_low/g;
tau_low=(R0_Full/R0_HH_low)*tau; % This is "naive" estimate which should be lower than one we want
tau_high=10*tau_low; % This should be high enough

flag=0;
while flag==0
    tau_now=0.5*(tau_high+tau_low);
    beta_ext=tau_now*d_int;
    beta_int=tau_now*d_ext;
    [Rs_HH_now, r_HH_now] = Grab_Epi_Data(E,E2,0,beta_ext,0, 0, beta_int, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
    if abs(r_HH_now-r_Full)<1e-4
        flag=1;
    end
    if r_HH_now<r_Full
        tau_low=tau_now;
    else
        tau_high=tau_now;
    end
end

tau_HH=tau_now;

Rs_HH=Rs_HH_now; r_HH=r_HH_now;

R0_HH=1+r_HH/g;

% Now do equilibrium prevalence
iter = 1; % Helps us keep track of how far along we are, if needed
P_R=0;
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_,P_R,H_Eq, nVect, nTicker, Q_int, Q_ext, Q_demo]=HH_demo_structured(1e-6,1e-6, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
        flag = 1;
    catch
        warning('error finding initial guess state')
    end
end

PQ_Sum = H_*E; % Demographic quantity needed for external infection probabilities
% Now do recursive process to find equilibrium
flag = 0;
while flag==0
    % Note that we have no structured infection here
    Unstr_Inf = I_bar*beta_ext*ones(1,maxN-1);
    try
        [I_bar_new,I_T_new,~,P_R_new,H_Eq_new, ~, ~, ~, ~, ~]=HH_demo_structured(Unstr_Inf,0, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
        if sum(abs(H_Eq-H_Eq_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_T_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-P_R_new)/P_R_new<1e-6
            flag=1;
        end
        I_bar = I_bar_new; I_T=I_T_new; P_R = P_R_new; H_Eq = H_Eq_new; % Get new ambient conditions
    catch
        warning(['error on iteration' num2str(iter)]);
    end
    iter = iter+1;
    if iter>100
        flag = 1;
        warning('Failed to converge for i=1');
    end
end
Prev_HH=I_bar; PR_HH=P_R;
Equil_HH=H_Eq;

filename=['ModelOutput/UKStructureEqDists_' datestr(now,'ddmmyy_HHMMSS')];
save(filename,'Equil_Full','Equil_HH','Equil_POLY','Equil_Hom');

filename=['ModelOutput/UKStructureOutputs_' datestr(now,'ddmmyy_HHMMSS')];
save(filename,'PR_Full','PR_Hom','PR_POLY','PR_HH',...
    'Prev_Full','Prev_Hom','Prev_POLY','Prev_HH',...
    'R0_Full','R0_Hom','R0_POLY','R0_HH',...
    'Rs_Full','Rs_Hom','Rs_POLY','Rs_HH',...
    'r_Full','r_Hom','r_POLY','r_HH',...
    'tau_Hom','tau_POLY','tau_HH');
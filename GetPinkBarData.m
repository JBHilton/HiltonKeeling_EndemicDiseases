% Calculates endemic equilibria for homogeneous mixing model, matching
% early growth rate with full model results

currenttime=datestr(now,'ddmmyy_HHMMSS');

%% First get UK demog data
load('Parameters/UK_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/UKDemography.mat')

d_int=sum(ClassProb*(D_All-D_Ext));
d_ext=sum(ClassProb*D_Ext);

%% UK, measles
g=1/7; % Infectious period
tau=-log(0.244)*g/d_int_kids; % Unit time transmission rate based on Hope-Simpson
beta_int=tau*d_int; % Internal contact rate
beta_ext=tau*d_ext; % Unstructured external contact rate

% Do scaling of WAIFW matrix by multiplying through by transmission rate
Inf_Ext=tau*D_Ext; % This is for HH-structured model where we just do external as polymod

% Calculate Rs and r so we can quote it (although we should have this in
% other results
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

beta_Hom=g*R0_Full; % Unstructured external contact rate

% Calculate Rs and r so we can quote it (although we should have this in
% other results
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
P_R=0; % We assume new parents have been subject to the vaccine program and form partnerships randomly
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
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

save(['Joes_Output/NumericalOutputs/UKMeaslesHomEqDist_' currenttime],'Equil_Hom');

%% Now get Kenya demog data
load('Parameters/Kenya_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/KenyaDemography.mat')


d_int=sum(ClassProb*(D_All-D_Ext));
d_ext=sum(ClassProb*D_Ext);

%% Kenya, measles
g=1/7; % Infectious period
tau=-log(0.244)*g/d_int_kids; % Unit time transmission rate based on Hope-Simpson
beta_int=tau*d_int; % Internal contact rate
beta_ext=tau*d_ext; % Unstructured external contact rate

% Do scaling of WAIFW matrix by multiplying through by transmission rate
Inf_Ext=tau*D_Ext; % This is for HH-structured model where we just do external as polymod

% Calculate Rs and r so we can quote it (although we should have this in
% other results
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

beta_Hom=g*R0_Full; % Unstructured external contact rate

% Calculate Rs and r so we can quote it (although we should have this in
% other results
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
P_R=0; % We assume new parents have been subject to the vaccine program and form partnerships randomly
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
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

save(['ModelOutput/KenyaMeaslesHomEqDist_' currenttime],'Equil_Hom');

%% Kenya, mumps
g=1/8; % Infectious period
tau=-log(0.689)*g/d_int_kids; % Unit time transmission rate based on Hope-Simpson
beta_int=tau*d_int; % Internal contact rate
beta_ext=tau*d_ext; % Unstructured external contact rate

% Do scaling of WAIFW matrix by multiplying through by transmission rate
Inf_Ext=tau*D_Ext; % This is for HH-structured model where we just do external as polymod

% Calculate Rs and r so we can quote it (although we should have this in
% other results
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

beta_Hom=g*R0_Full; % Unstructured external contact rate

% Calculate Rs and r so we can quote it (although we should have this in
% other results
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
P_R=0; % We assume new parents have been subject to the vaccine program and form partnerships randomly
flag = 0;
while flag==0
    try
        [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
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

save(['ModelOutput/KenyaMumpsHomEqDist_' currenttime],'Equil_Hom');
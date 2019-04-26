% Does the following calculations for a disease in the UK with g=1/7 for
% all four transmission structures:
%   1: r as a function of tau
%   2: R_* as a function of tau
%   3: Endemic prevalence as a function of tau
%   4: PR as a function of tau

load('Parameters/UK_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/UKDemography.mat')

% Define infectious parameters
g=1/7;
d_int=sum(ClassProb*(D_All-D_Ext));
d_ext=sum(ClassProb*D_Ext);

% Cycle over sigma values and calculate early growth parameters,
% prevalence, prob of acquiring during childhood.
tau=[0.4,1.8]; % The results for sigma=0:0.1:1 are plotted in the paper - this takes a long time to calculate for the Kenya-like population so we just do a mumps-like and a measles-like value for this example
Rs_Hom=zeros(1,length(tau)); r_Hom=zeros(1,length(tau));
Prev_Hom=zeros(1,length(tau)); PR_Hom=zeros(1,length(tau));
Rs_POLY=zeros(1,length(tau)); r_POLY=zeros(1,length(tau));
Prev_POLY=zeros(1,length(tau)); PR_POLY=zeros(1,length(tau));
Rs_HH=zeros(1,length(tau)); r_HH=zeros(1,length(tau));
Prev_HH=zeros(1,length(tau)); PR_HH=zeros(1,length(tau));
Rs_Full=zeros(1,length(tau)); r_Full=zeros(1,length(tau));
Prev_Full=zeros(1,length(tau)); PR_Full=zeros(1,length(tau));

for i=1:length(tau)
    beta_int=tau(i)*d_int; % Internal contact rate
    beta_ext=tau(i)*d_ext; % Unstructured external contact rate
    beta_all=tau(i)*(d_int+d_ext); % Non-age-structured contact rate
    
    % Do scaling of WAIFW matrix by multiplying through by transmission rate
    Inf_Ext=tau(i)*D_Ext; % This is for HH-structured model where we just do external as polymod
    Inf_All=tau(i)*D_All; % This is for model without HH structure where we just do age structure
    % First calculate early growth parameters
    
    % Homogeneous model
    starttime=cputime;
    flag = 0;
    while flag==0
        try
            [Rs_Hom(i), r_Hom(i)] = Grab_Epi_Data(E,E2,Inf_Ext,beta_all,0, 0, 0, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
            flag = 1;
        catch
            warning('Error in Grab_Epi_Data.m');
        end
    end
    
    % Purely age-structured model
    flag = 0;
    while flag==0
        try
            [Rs_POLY(i), r_POLY(i)] = Grab_Epi_Data(E,E2,Inf_All,0,1, 0, 0, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
            flag = 1;
        catch
            warning('Error in Grab_Epi_Data.m');
        end
    end
    
    % Purely household-structured model
    flag = 0;
    while flag==0
        try
            [Rs_HH(i), r_HH(i)] = Grab_Epi_Data(E,E2,Inf_Ext,beta_ext,0, 0, beta_int, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
            flag = 1;
        catch
            warning('Error in Grab_Epi_Data.m');
        end
    end
    
    % Full model
    flag = 0;
    while flag==0
        try
            [Rs_Full(i), r_Full(i)] = Grab_Epi_Data(E,E2,Inf_Ext,beta_ext,1, 0, beta_int, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
            flag = 1;
        catch
            warning('Error in Grab_Epi_Data.m');
        end
    end
    
    % Now do equilibrium prevalence
    
    % Homogeneous model
    if r_Hom(i)>0
        iter = 1; % Helps us keep track of how far along we are, if needed
        P_R=0; 
        flag = 0;
        while flag==0
            try
                [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
                flag = 1;
            catch
                warning('error finding initial guess state')
            end
        end
        
        PQ_Sum = H_T*E; % Demographic quantity needed for external infection probabilities
        % Now do recursive process to find equilibrium
        flag = 0;
        while flag==0
            Unstr_Inf = I_bar*beta_all*ones(1,maxN-1); % Unstructured ext. inf.
            try
                [I_bar_new,I_ave_new,~,R_leave_new,FD_new, ~, ~, ~, ~, ~]=HH_demo_structured(Unstr_Inf,0, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
                if sum(abs(H_Eq-FD_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_ave_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-R_leave_new)/R_leave_new<1e-6
                    flag=1;
                end
                I_bar = I_bar_new; I_T=I_ave_new; P_R = R_leave_new; H_Eq = FD_new; % Get new ambient conditions
            catch
                warning(['error on iteration' num2str(iter)]);
            end
            iter = iter+1;
            if iter>100
                flag = 1;
                warning('Failed to converge for i=1');
            end
        end
        Prev_Hom(i)=I_bar; PR_Hom(i)=P_R;
    else
        Prev_Hom(i)=0; PR_Hom(i)=0;
    end
    
    % Purely age-structured model
    if r_POLY(i)>0
        iter = 1; % Helps us keep track of how far along we are, if needed
        P_R=0; 
        flag = 0;
        while flag==0
            try
                [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, beta_all, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
                flag = 1;
            catch
                warning('error finding initial guess state')
            end
        end
        
        PQ_Sum = H_T*E; % Demographic quantity needed for external infection probabilities
        % Now do recursive process to find equilibrium
        flag = 0;
        while flag==0
            Struc_Inf = E*Inf_All*E2*I_T'; % Age-structured ext. inf.
            try
                [I_bar_new,I_ave_new,~,R_leave_new,FD_new, ~, ~, ~, ~, ~]=HH_demo_structured(0,Struc_Inf, P_R, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
                if sum(abs(H_Eq-FD_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_ave_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-R_leave_new)/R_leave_new<1e-6
                    flag=1;
                end
                I_bar = I_bar_new; I_T=I_ave_new; P_R = R_leave_new; H_Eq = FD_new; % Get new ambient conditions
            catch
                warning(['error on iteration' num2str(iter)]);
            end
            iter = iter+1;
            if iter>100
                flag = 1;
                warning('Failed to converge for i=1');
            end
        end
        Prev_POLY(i)=I_bar; PR_POLY(i)=P_R;
    else
        Prev_POLY(i)=0; PR_POLY(i)=0;
    end
    
    % Purely household-structured model
    if r_HH(i)>0
        iter = 1; % Helps us keep track of how far along we are, if needed
        P_R=0; 
        flag = 0;
        while flag==0
            try
                [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
                flag = 1;
            catch
                warning('error finding initial guess state')
            end
        end
        
        PQ_Sum = H_T*E; % Demographic quantity needed for external infection probabilities
        % Now do recursive process to find equilibrium
        flag = 0;
        while flag==0
            Unstr_Inf = I_bar*beta_ext*ones(1,maxN-1);
            try
                [I_bar_new,I_ave_new,~,R_leave_new,FD_new, ~, ~, ~, ~, ~]=HH_demo_structured(Unstr_Inf,0, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
                if sum(abs(H_Eq-FD_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_ave_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-R_leave_new)/R_leave_new<1e-6
                    flag=1;
                end
                I_bar = I_bar_new; I_T=I_ave_new; P_R = R_leave_new; H_Eq = FD_new; % Get new ambient conditions
            catch
                warning(['error on iteration' num2str(iter)]);
            end
            iter = iter+1;
            if iter>100
                flag = 1;
                warning('Failed to converge for i=1');
            end
        end
        Prev_HH(i)=I_bar; PR_HH(i)=P_R;
    else
        Prev_HH(i)=0; PR_HH(i)=0;
    end
    
    % Full model
    if r_Full(i)>0
        iter = 1; % Helps us keep track of how far along we are, if needed
        P_R=0; 
        flag = 0;
        while flag==0
            try
                [I_bar,I_T,H_T,P_R,H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(1e-6,1e-6, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,DiseaseFree);
                flag = 1;
            catch
                warning('error finding initial guess state')
            end
        end
        
        PQ_Sum = H_T*E; % Demographic quantity needed for external infection probabilities
        % Now do recursive process to find equilibrium
        flag = 0;
        while flag==0
            Struc_Inf = E*Inf_Ext*E2*I_T'; % Age-structured ext. inf.
            try
                [I_bar_new,I_ave_new,~,R_leave_new,FD_new, ~, ~, ~, ~, ~]=HH_demo_structured(0,Struc_Inf, P_R, beta_int, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, H_Eq);
                if sum(abs(H_Eq-FD_new))<length(H_Eq)*1e-6 & sum(abs(I_T-I_ave_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (P_R-R_leave_new)/R_leave_new<1e-6
                    flag=1;
                end
                I_bar = I_bar_new; I_T=I_ave_new; P_R = R_leave_new; H_Eq = FD_new; % Get new ambient conditions
            catch
                warning(['error on iteration' num2str(iter)]);
            end
            iter = iter+1;
            if iter>100
                flag = 1;
                warning('Failed to converge for i=1');
            end
        end
        Prev_Full(i)=I_bar; PR_Full(i)=P_R;
    else
        Prev_Full(i)=0; PR_Full(i)=0;
    end
    
    fintime=cputime-starttime;
    disp(['Completed calculation ' num2str(i) ' of ' num2str(length(tau)) ' in ' num2str(fintime) ' seconds.']);
end
filename=['ModelOutput/UKTau7Loop_' datestr(now,'ddmmyy_HHMMSS')];
save(filename, 'tau', 'PR_Hom', 'PR_POLY', 'PR_HH', 'PR_Full',...
    'Prev_Hom', 'Prev_POLY', 'Prev_HH', 'Prev_Full',...
    'r_Hom', 'r_POLY', 'r_HH', 'r_Full',...
   'Rs_Hom', 'Rs_POLY', 'Rs_HH', 'Rs_Full');
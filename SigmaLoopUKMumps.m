% Does the following calculations for mumps in the UK:
%   1: r as a function of sigma
%   2: R_* as a function of sigma
%   3: Endemic prevalence as a function of sigma
%   4: PR as a function of sigma

load('Parameters/UK_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/UKDemography.mat')

% Define infectious parameters
g=1/8;
d_int=sum(ClassProb*(D_All-D_Ext));
d_ext=sum(ClassProb*D_Ext);
tau=-log(0.689)*g/Ave_int_dur_kids; % Unit time transmission rate based on Hope-Simpson
save('Parameters/MumpsTransRate.mat','tau'); % Save this to use in Kenya calculations
beta_int=tau*d_int; % Internal contact rate
beta_ext=tau*d_ext; % Unstructured external contact rate
beta_all=tau*(d_int+d_ext); % Non-age-structured contact rate

% Do scaling of contact matrix by multiplying through by transmission rate
Inf_Ext=tau*D_Ext; % This is for HH-structured model where we just do external as polymod
Inf_All=tau*D_All; % This is for model without HH structure where we just do age structure

% Cycle over sigma values and calculate early growth parameters,
% prevalence, prob of acquiring during childhood.
sigma=0:0.1:1.5;
Rs_Full=zeros(1,length(sigma)); r_Full=zeros(1,length(sigma));
Prev_Full=zeros(1,length(sigma)); PR_Full=zeros(1,length(sigma));
Rs_AgeOnly=zeros(1,length(sigma)); r_AgeOnly=zeros(1,length(sigma));
Prev_AgeOnly=zeros(1,length(sigma)); PR_AgeOnly=zeros(1,length(sigma));

for i=1:length(sigma)
    % First calculate early growth parameters
    % Full model
    starttime=cputime;
    flag = 0;
    while flag==0
        try
            [Rs_Full(i), r_Full(i)] = Grab_Epi_Data(Q,Q2,Inf_Ext,beta_ext,sigma(i), 0, beta_int, g, kB,kL,kB,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
            flag = 1;
        catch
            warning('Error in Grab_Epi_Data.m for unvaccinated population');
        end
    end
    
    % Age structure only
    flag = 0;
    while flag==0
        try
            [Rs_AgeOnly(i), r_AgeOnly(i)] = Grab_Epi_Data(Q,Q2,Inf_All,beta_all,sigma(i), 0, 0, g, kB,kL,kB,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid);
            flag = 1;
        catch
            warning('Error in Grab_Epi_Data.m for unvaccinated population');
        end
    end
    
    % Now do equilibrium prevalence
    % Full model
    if r_Full(i)>0
        iter = 1; % Helps us keep track of how far along we are, if needed
        R_leave=0; % We assume new parents have been subject to the vaccine program and form partnerships randomly
        flag = 0;
        while flag==0
            try
                [I_bar,I_ave,~,R_leave,FD, nVect, nTicker, N, nM, eM]=HH_demo_structured(1e-6,1e-6, R_leave, beta_int, g, 0, kB,kL,kB,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0,DiseaseFree);
                flag = 1;
            catch
                warning('error finding initial guess state')
            end
        end
        
        % Now do recursive process to find equilibrium
        flag = 0;
        while flag==0
            Unstr_Inf = I_bar*beta_ext*sign(1-sigma(i))*(1-sigma(i))*ones(1,maxN-1); % Unstructured ext. inf.
            Struc_Inf = min(sigma(i),1)*Q*Inf_Ext*Q2*I_ave'; % Age-structured ext. inf.
            try
                [I_bar_new,I_ave_new,~,R_leave_new,FD_new, ~, ~, ~, ~, ~]=HH_demo_structured(Unstr_Inf,Struc_Inf, R_leave, beta_int, g, 0, kB,kL,kB,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0, 0, FD);
                if sum(abs(FD-FD_new))<length(FD)*1e-6 & sum(abs(I_ave-I_ave_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (R_leave-R_leave_new)/R_leave_new<1e-6
                    flag=1;
                end
                I_bar = I_bar_new; I_ave=I_ave_new; R_leave = R_leave_new; FD = FD_new; % Get new ambient conditions
            catch
                warning(['error on iteration' num2str(iter)]);
            end
            iter = iter+1;
            if iter>100
                flag = 1;
                warning('Failed to converge for i=1');
            end
        end
        Prev_Full(i)=I_bar; PR_Full(i)=R_leave;
    else
        Prev_Full(i)=0; PR_Full(i)=0;
    end
    if abs(sigma(i)-1)<0.01
        Equil_MumpsUK=FD;
        filename=['NumericalOutputs/UKMumpsEqDist_' datestr(now,'ddmmyy_HHMMSS')];
        save(filename,'Equil_MumpsUK');
    end
    
    % Age only
    if r_AgeOnly(i)>0
        iter = 1; % Helps us keep track of how far along we are, if needed
        R_leave=0; % We assume new parents have been subject to the vaccine program and form partnerships randomly
        flag = 0;
        while flag==0
            try
                [I_bar,I_ave,~,R_leave,FD, nVect, nTicker, N, nM, eM]=HH_demo_structured(1e-6,1e-6, R_leave, beta_all, g, 0, kB,kL,kB,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0,DiseaseFree);
                flag = 1;
            catch
                warning('error finding initial guess state')
            end
        end
        
        % Now do recursive process to find equilibrium
        flag = 0;
        while flag==0
            Unstr_Inf = I_bar*beta_all*sign(1-sigma(i))*(1-sigma(i))*ones(1,maxN-1); % Unstructured ext. inf.
            Struc_Inf = min(sigma(i),1)*Q*Inf_All*Q2*I_ave'; % Age-structured ext. inf.
            try
                [I_bar_new,I_ave_new,~,R_leave_new,FD_new, ~, ~, ~, ~, ~]=HH_demo_structured(Unstr_Inf,Struc_Inf, R_leave, 0, g, 0, kB,kL,kB,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0, FD);
                if sum(abs(FD-FD_new))<length(FD)*1e-6 & sum(abs(I_ave-I_ave_new))<1e-6 & (I_bar-I_bar_new)/I_bar_new<1e-6 & (R_leave-R_leave_new)/R_leave_new<1e-6
                    flag=1;
                end
                I_bar = I_bar_new; I_ave=I_ave_new; R_leave = R_leave_new; FD = FD_new; % Get new ambient conditions
            catch
                warning(['error on iteration' num2str(iter)]);
            end
            iter = iter+1;
            if iter>100
                flag = 1;
                warning('Failed to converge for i=1');
            end
        end
        Prev_AgeOnly(i)=I_bar; PR_AgeOnly(i)=R_leave;
    else
        Prev_AgeOnly(i)=0; PR_AgeOnly(i)=0;
    end
    fintime=cputime-starttime;
    disp(['Completed calculation ' num2str(i) ' of ' num2str(length(sigma)) ' in ' num2str(fintime) ' seconds.']);
end
filename=['NumericalOutputs/UKMumpsSigmaLoop_' datestr(now,'ddmmyy_HHMMSS')];
save(filename, 'sigma', 'PR_AgeOnly', 'PR_Full', 'Prev_AgeOnly', 'Prev_Full', 'r_AgeOnly', 'r_Full', 'Rs_AgeOnly', 'Rs_Full');
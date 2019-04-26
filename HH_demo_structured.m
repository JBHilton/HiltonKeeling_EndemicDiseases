function [I_bar, I_T, H_T, P_R, H_Eq, nVect, nTicker, Q_demo, Q_int, Q_ext]=HH_demo_structured(Unstr_Inf,Struc_Inf, P_R, beta, g, Vacc_at_birth, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid,DISP,eig_0)
% HH_demo_structured calculates the transition matrix Q=(Q_demo+Q_int+Q_ext) of the
% household model with demography along with its equilibrium distribution
% and important statistics of this distribution. Its inputs are the
% parameters of the model along with some index vectors and an ansatz for
% the eigenvalue calculation.
%   All parameters and index vectors are defined in the scripts that call
%   this one. The ansatz eig_0 is intended to be an already known
%   equilibrium distribution for a set of parameters near to the ones we
%   are currently using - i.e. the previous equilibrium if we are changing
%   parameters iteratively. If no such distribution is known, set eig_0=0.
%   The DISP flag decides whether to display or supress details about the
%   eigenvalue calculation. Set it to zero to supress, otherwise it will
%   display the tolerance and starting ansatz used.

maxN = find(StopProb==1,1);
if isempty(maxN)
    StopProb(end)=1; % Make sure stopping probabilities are well-defined!
    maxN=length(StopProb);
end

%% SIR classes without external transmission.

[S,I,R]=meshgrid([0:maxN],[0:maxN],[0:maxN]); m=find(S+I+R<=maxN & S+I+R>=2)'; % Impose household size constraint
VectN=S(m)+I(m)+R(m); VectS=S(m); VectI=I(m); VectR=R(m); Vect=[VectS; VectI; VectR]; % Ordering of all states
Q_int_short=sparse(1,1,0,length(VectN),length(VectN)); % Transition matrix for internal epidemic
for LeaveStates=1:length(VectN)
    R=beta*VectS(LeaveStates)*VectI(LeaveStates)/(VectN(LeaveStates)-1) + g*VectI(LeaveStates); Q_int_short(LeaveStates,LeaveStates)=-R;
    R=beta*VectS(LeaveStates)*VectI(LeaveStates)/(VectN(LeaveStates)-1); % Infection
    if R>0
        q=find(VectS==VectS(LeaveStates)-1 & VectI==VectI(LeaveStates)+1 & VectR==VectR(LeaveStates));
        Q_int_short(LeaveStates,q)=Q_int_short(LeaveStates,q)+R; % Add to transition matrix
    end
    R=g*VectI(LeaveStates); % Recovery
    if R>0
        q=find(VectS==VectS(LeaveStates) & VectI==VectI(LeaveStates)-1 & VectR==VectR(LeaveStates)+1);
        Q_int_short(LeaveStates,q)=Q_int_short(LeaveStates,q)+R;
    end
end

%%
% The Prob_at_(x,y) terms are the probability of going to ticker location x
% given that the household is size y.

%Set up movement rates through waiting for birth classes
% k is shape parameter (number of stages) in gamma waiting process
Birth_Ticker_Rates=ones(maxN,kB)*kB/TB;
%Set up what to do after a birth
% Prob_at_Birth is (1+kB)x(maxN) determines the state you go to
% after a birth -- either back to ticker=0 or on to ticker=1+kB
Prob_at_Birth=sparse(1+kB+0*StopProb, 1:maxN, StopProb, 1+kB, maxN);
Prob_at_Birth=Prob_at_Birth+sparse(1+0*StopProb, 1:maxN, 1-StopProb, 1+kB,maxN);
Prob_at_Birth(1,maxN)=0; Prob_at_Birth(1+kB,maxN)=1;


%Set up movement rates through waiting to leave home classes
%note that this is for the FIRST child to leave; subsequent children follow
%the birth timings.
% Expected first leave wait time is (mean age at leaving - (number of
% younger siblings)*(time between births)), i.e. counts from last birth.
% This could be negative - this makes the script return error
Leave_Ticker_Rates=ones(maxN,kL)*kL./(TL - ([1:maxN]'-3)*ones(1,kL)*TB);
Prob_at_Leaving=sparse(1+kB+kL+0*[4:maxN],[4:maxN],1+0*[4:maxN],1+kB+kL+kB,maxN);
Prob_at_Leaving=Prob_at_Leaving+sparse(1+kB+kL+kB,3,1,1+kB+kL+kB,maxN);

% Expected intervals for subsequent leavings are just birth intervals
Leave_Ticker_Rates2=ones(maxN,kB)*kB./(TB);
Leave_Ticker_Rates2(end,:)=0; % as cannot be in state maxN as one child has left
%Set up what to do after a birth
Prob_at_Leaving2=sparse(1+kB+kL+0*[4:maxN],[4:maxN],1+0*[4:maxN],1+kB+kL+kB,maxN);
Prob_at_Leaving2=Prob_at_Leaving2+sparse(1+kB+kL+kB,3,1,1+kB+kL+kB,maxN);

%Set up movement rates through waiting to be re-born classes
Reset_Ticker_Rates=ones(maxN,kR)*kR./(TR);
%Set up what to do after a reset
if max(VectR)<2 % Inf_at_Reset is probability of a child having been infected before starting new home
    Inf_at_Reset=0;
else
    Inf_at_Reset=P_R; % Proportion of newly recovered households
end
Ticker_at_Reset(1)=1-StopProb(2); Ticker_at_Reset(1+kB+kL+kB)=StopProb(2);
% Instant reset if no kids, otherwise reset at end
%%

% Add_Demography calculates Q_demo and repeats all the infectious status
% vectors/transition matrix so that it covers all the demographic classes.

[nVect,nVectS,nVectI,nVectR,nVectN,nTicker,Q_int,Q_demo]=Add_Demography(Vect,VectS,VectI,VectR,VectN,Q_int_short, Vacc_at_birth, Birth_Ticker_Rates,Prob_at_Birth, Leave_Ticker_Rates, Prob_at_Leaving, Leave_Ticker_Rates2, Prob_at_Leaving2, Reset_Ticker_Rates, Inf_at_Reset, Ticker_at_Reset);

% We now calculate eM, the transition matrix of external infection events.

 T = nVectN*(kB+kL+kB+kR)+nTicker; % Full demo state of each demo-epi state

if size(Struc_Inf)==1
    Struc_Inf=Struc_Inf*ones(1,length(DemGrid));
end
if size(Unstr_Inf)==1
    Unstr_Inf=Unstr_Inf*ones(1,maxN-1);
end

Q_ext=sparse(1,1,0,length(nVectN),length(nVectN)); % eM is external infection
for LeaveStates=1:length(nVectN)
    R=Unstr_Inf(nVectN(LeaveStates)-1)+Struc_Inf(find(DemGrid==T(LeaveStates)))/nVectN(LeaveStates);
    R=R*nVectS(LeaveStates); %%% This assumes external is per person %%%
    if R>0
        Q_ext(LeaveStates,LeaveStates)=-R;
        q=find(nVectS==nVectS(LeaveStates)-1 & nVectI==nVectI(LeaveStates)+1 & nVectR==nVectR(LeaveStates) & nTicker==nTicker(LeaveStates));
        Q_ext(LeaveStates,q)=Q_ext(LeaveStates,q)+R;
    end
end


%%
% Now calculate the zero-eigenvector (largest right eigenvector) to get the
% equilibrium distribution. If no eig_0 is given, we just calculate the
% disease-free equilibrium (using two adults with no children as starting
% guess) and add a small amount of infection to give us our starting point.
% If this fails, we use a random starting point. If eig_0 is given, we try
% this first.



userconv=0; % This flag tells us if the user-defined ansatz led to convergence
if length(eig_0)==length(nVectS) % Check if we have valid ansatz
    useropts.v0=eig_0; % This is the user-defined ansatz
    useropts.tol=1e-10; % This is Matlab's default tolerance - we will lower it if it doesn't work
    tolflag=0;
    while (userconv==0)&&(tolflag==0)
        try
            [H_Eq,~]=eigs((Q_int+Q_demo+Q_ext)',1,'lr',useropts);
            userconv=1;
            if DISP>0
                disp(['Found equilibrium with user ansatz to tolerance ' num2str(useropts.tol)]);
            end
        catch            
            useropts.tol=10*useropts.tol; % Increase tolerance if we don't converge
            if useropts.tol>1e-3 % Give up if tolerance gets too high
                tolflag=1;
                if DISP>0
                    warning('User-defined ansatz failed, moving on to default ansatz');
                end
            end
        end
    end
end

defconv=0; % This flag tells us if default ansatz converges
if userconv==0 % Do this if we didn't get an eigenvector from user input
    % Define default starting conditions for system by doing demographic
    % equilibrium with a bit of infection
    demodist=zeros(length(nVectN),1);
    demodist(nTicker==1&nVectS==2)=(1-P_R)^2;
    demodist(nTicker==1&nVectS==1&nVectR==1)=P_R*(1-P_R);
    demodist(nTicker==1&nVectR==2)=P_R^2;
    opts.v0 = demodist;
    [demodist,~]=eigs(Q_demo',1,'lr',opts);
    demodist=demodist/sum(demodist);
    demodist(demodist<0)=0;
    demodist(nVectI==1)=1e-6;
    demodist=demodist/sum(demodist);
    defopts.v0=demodist; % This is our default ansatz
    defopts.tol=1e-10; % This is Matlab's default tolerance
    tolflag=0;
    while (defconv==0)&&(tolflag==0)
        try
            [H_Eq,~]=eigs((Q_int+Q_demo+Q_ext)',1,'lr',defopts);
            defconv=1;
            if DISP>0
                disp(['Found equilibrium with default ansatz to tolerance ' num2str(defopts.tol)]);
            end
        catch
            defopts.tol=10*defopts.tol; % Increase tolerance if we don't converge
            if defopts.tol>1e-3 % Give up if tolerance gets too high
                tolflag=1;
                if DISP>0
                    warning('Default ansatz failed, moving on to random ansatz');
                end;
            end
        end
    end
end

randconv=0; % If user-defined and default ansatzes don't work, we do random
if (userconv==0)&&(defconv==0)
    randopts.tol=1e-6; % Start from a bigger tolerance since we already know it's hard
    tolflag=0;
    while (randconv==0)&&(tolflag==0) % No tolflag, we just keep guessing in this case
        try
            [H_Eq,~]=eigs((Q_int+Q_demo+Q_ext)',1,'lr',randopts);
            randconv=1;
            if DISP>0
                disp(['Found equilibrium with random ansatz to tolerance ' num2str(randopts.tol)]);
            end
        catch
            randopts.tol=10*randopts.tol; % Increase tolerance if we don't converge
            if randopts.tol<1e-3 % Once we get to 1e-3 just keep trying over and over
                tolflag=1;
                if DISP>0
                    warning('Random ansatz, change to smallest magnitude option');
                end
            end
        end
    end
end

smconv=0; % If user-defined and default ansatzes don't work, we do random
if (userconv==0)&&(defconv==0)
    smopts.tol=1e-6; % Start from a bigger tolerance since we already know it's hard
    while (smconv==0) % No tolflag, we just keep guessing in this case
        try
            [H_Eq,~]=eigs((Q_int+Q_demo+Q_ext)',1,'sm',smopts);
            smconv=1;
            if DISP>0
                disp(['Found equilibrium with random ansatz and smallest magnitude option to tolerance ' num2str(smopts.tol)]);
            end
        catch
            if DISP>0
                warning('Random ansatzand smallest magnitude option failed to converge');
            end
            if smopts.tol<1e-3 % Once we get to 1e-3 just keep trying over and over
                smopts.tol=10*smopts.tol; % Increase tolerance if we don't converge
            end
        end
    end
end

H_Eq=H_Eq/sum(H_Eq); % Normalise

I_bar=sum(nVectI.*H_Eq')/sum(nVectN.*H_Eq'); % Population-level infectious prevalence

% Next loop calculates demographic class-stratified infectious prevalence
H_T = 0*NGrid;
Index=cell(1,length(NGrid));
I_T=zeros(1,length(NGrid));
for i=1:length(NGrid)
    Index{i}=find(nVectN==NGrid(i)&nTicker==tickGrid(i));
    H_T(i)=sum(H_Eq(Index{i}));
    I_T(i)=nVectI(Index{i})*H_Eq(Index{i})/(H_T(i)*NGrid(i));
end

% Finally, calculate P_R
LeaveStates=find(nTicker==(kB+kL) | nTicker==(kB+kL+kB)); % States at leaving stage
nER_Start=2*P_R;
nER_Leave = nVectR; % Total exposure level at leaving
P_R = sum( ((nER_Leave(LeaveStates) - nER_Start)./(nVectN(LeaveStates)-nER_Start)).*H_Eq(LeaveStates)')/sum(H_Eq(LeaveStates));
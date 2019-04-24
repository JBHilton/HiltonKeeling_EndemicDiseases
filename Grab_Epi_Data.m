function [Rstar, r] = Grab_Epi_Data(E,E2,D, alpha, sigma, P_R, beta, g, kB,kL,kR,TB,TL,TR,StopProb,DemGrid,NGrid,tickGrid)
% Grab_Epi_Data calculates the household-level reproductive number Rstar
% and the early exponential growth rate r for the household model with
% demography in the absence of vaccination.
%   The meaning of the inputs is explained in the scripts that call
%   this function.


maxN = find(StopProb==1,1);
if isempty(maxN)
    StopProb(end)=1; % Make sure stopping probabilities are well-defined!
    maxN=length(StopProb);
end

% Now calculate disease-free equilibrium distribution and internal
% transmission matrix
[~, ~, H_T, ~, H_Eq, nVect, nTicker, ~,~, ~]=HH_demo_structured(0, 0, 0, 0, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,0);
[~, ~, ~, ~, ~, ~, ~, ~, Q_int]=HH_demo_structured(0, 0, P_R, beta, g, 0, kB,kL,kR,TB,TL,TR,StopProb,NGrid,tickGrid,DemGrid, 0,H_Eq);
Infected=find(nVect(2,:)+nVect(3,:)>0); % Indices of states with infection history
nVectN=sum(nVect,1); PopSize=nVectN*H_Eq; %  PopSize is mean household size

T = nVectN*(kB+kL+kB+kR)+nTicker; % Indices of full demographic states

F_D=H_Eq; F_D(F_D<0)=0; F_D(nVect(2,:)==1&nVect(3,:)==0)=1e-5;  % add a bit of infection

Q_int_long=Q_int;
m=find(nVect(2,:)+nVect(3,:)==0);
for i=1:length(m)
    Q_int_long(m(i),m(i))=0; % Get rid of outflow so that we can study growth in infected states
end

% In this section we calculate exponential growth rate r=\dot{I}/I.
% We repeatedly run the linearised equations forward over short distances,
% scaling down the additional infection at the end of each run. In this way
% we get the quasi-equilibrium distribution of infection in the early
% stages of the epidemic, which is the eigenvector associated with the
% eigenvalue r, which is the growth rate of the system.
flag=1;
attempts=0;
while flag
    [t, X] = ode45( @(t,p) ODEs(t,p,nVect, nTicker, Q_int_long, DemGrid, T, H_T, E, E2, alpha, D, sigma, kB+kL+kB+kR, PopSize) , [0 50], F_D, odeset('RelTol',1e-8, 'AbsTol', 1e-10,'NonNegative',[1:length(F_D)]) );
    FFD=X(end,:);
    FFD(nVect(2,:)+nVect(3,:)>0)=...
        1e-5*FFD(nVect(2,:)+nVect(3,:)>0)/sum(FFD(nVect(2,:)+nVect(3,:)>0));
    dX=ODEs(0,X(end,:)',nVect, nTicker, Q_int_long, DemGrid, T, H_T, E, E2, alpha, D, sigma, kB+kL+kB+kR, PopSize);
    Xnow=X(end,:)';
    R=dX./Xnow; R=R((Xnow'>0)&(nVect(2,:)>0));
    attempts=attempts+1;
    sXn=sum(Xnow(nVect(2,:)>0)); sFD=sum(F_D(nVect(2,:)>0));
    fprintf(1,'Attempt %d) X=%g R=%g var(R)=%g [sum(end)=%g sum(start)=%g %g]\n',attempts,mean(Xnow),mean(R),var(R),sXn,sFD,log(sXn/sFD)/t(end));
    if var(R)<1e-7||attempts>99
        flag=0;
        if attempts>99
            R=log(sXn/sFD)/t(end);
        end
    end
    if sXn<sFD*5
        FFD(nVect(2,:)+nVect(3,:)>0)=1e-5*FFD(nVect(2,:)+nVect(3,:)>0)/sum(FFD(nVect(2,:)>0));
        if sXn<sFD*0.1 && attempts>1 && flag
            flag=0;
            R=log(sXn/sFD)/t(end);
        end
        if sXn<sFD*10 && attempts>10 && flag
            flag=0;
            R=log(sXn/sFD)/t(end);
        end
    end
    F_D=FFD;
end
r=mean(R); % Converged per capita growth rate

%%
% Now do R_star
%%%


Z=X(end,:)'/sum(X(end,Infected)); % Equilibrium distribution conditional on having infectious history

% Now get external infection matrix again

nVectI=nVect(2,:); nVectS=nVect(1,:);  nVectR=nVect(3,:);
nVectN=sum(nVect,1);
 
I_bar=sum(nVectI.*Z')/PopSize; % Population-level infectious prevalence
I_T=zeros(1,length(H_T));

for i=1:length(nVectN)
    State_i = find(DemGrid==T(i));
    I_T(State_i) = I_T(State_i) + (nVectI(i)/nVectN(i))*Z(i); % Demostate-stratified prevalence
end
I_T=I_T./H_T;

PQ_Sum = H_T*E;
Unstr_Inf = I_bar*alpha*sign(1-sigma)*(1-sigma)*ones(1,maxN-1); % Unstructured ext. inf.
Struc_Inf = min(sigma,1)*E*D*E2*I_T'; % Age-structured ext. inf.

Q_ext=sparse(1,1,0,length(nVectN),length(nVectN)); % Q_ext is external infection
for n=1:length(nVectN)
    if (nVect(2,n)+nVect(3,n))==0  %% as assuming infection is rare only to all S's - household level branching
        R=Unstr_Inf(nVectN(n)-1)+Struc_Inf(find(DemGrid==nVectN(n)*(kB+kL+kB+kR)+nTicker(n)))/nVectN(n);
        R=R*nVectS(n); %%% This assumes external is per person %%%
        if R>0
            q=find(nVectS==nVectS(n)-1 & nVectI==nVectI(n)+1 & nVectR==nVectR(n) & nTicker==nTicker(n));
            Q_ext(n,q)=Q_ext(n,q)+R;
        end
    end
end

% these are the new cases made by the eigen-distribution.
New_Cases=(Q_ext)'*H_Eq; % Rate of new case production
New_Cases=New_Cases/sum(New_Cases);

% Run the dynamics forward -- but NEW households are counted separately

% Transition matrix has extra dimension corresponding to "new house" state
Q_int_long(end+1,end+1)=0;
F_D=H_Eq+New_Cases; F_D(end+1)=0;
[t, X] = ode45( @(t,p) Extended_ODEs(t,p,nVect, nTicker, Q_int_long, DemGrid, T, H_T, E, E2, alpha, D, sigma, kB+kL+kB+kR, PopSize) , [0 50], F_D, odeset('RelTol',1e-8, 'AbsTol', 1e-10) );

Rstar=X(end,end); % Expected number of times naive household gets hit
end

function [dX] = ODEs(t, Z, nVect, nTicker, Q_int, DemGrid, T, H_T, E, E2, alpha, D, sigma, sumK, POP)
% This does the master equations except with diagonal set to zero, so
% probability mass accumulates in states rather than flowing out. This
% allows us to measure the growth in infection in the early stages of an
% outbreak.

nVectI=nVect(2,:); nVectS=nVect(1,:);  nVectR=nVect(3,:);
nVectN=sum(nVect,1); maxN=max(nVectN);
 
I_bar=sum(nVectI.*Z')/POP;
I_T=zeros(1,length(H_T));

for i=1:length(nVectN)
    State_i = find(DemGrid==T(i));
    I_T(State_i) = I_T(State_i) + (nVectI(i)/nVectN(i))*Z(i);
end
I_T=I_T./H_T;

PQ_Sum = H_T*E;
Unstr_Inf = I_bar*alpha*sign(1-sigma)*(1-sigma)*ones(1,maxN-1); % Unstructured ext. inf.
Struc_Inf = min(sigma,1)*E*D*E2*I_T'; % Age-structured ext. inf.

Q_ext=sparse(1,1,0,length(nVectN),length(nVectN));
for n=1:length(nVectN) % Builds matrix of external infection rates
    if (nVect(2,n)+nVect(3,n))==0  %% as assuming infection is rare only to all S's (branching assumption - Joe)
        R=Unstr_Inf(nVectN(n)-1)+Struc_Inf(find(DemGrid==nVectN(n)*(sumK)+nTicker(n)))/nVectN(n);
        R=R*nVectS(n); %%% This assumes external is per person %%%
        if R>0
            q=find(nVectS==nVectS(n)-1 & nVectI==nVectI(n)+1 & nVectR==nVectR(n) & nTicker==nTicker(n));
            Q_ext(n,q)=Q_ext(n,q)+R;
        end
    end
end

dX=(Q_int+Q_ext)'*Z;

end


% As above, but all newly infected households go to a new count -
% BRANCHING!
function [dX] = Extended_ODEs(t, ZZ, nVect, nTicker, Q_int, DemGrid, T, H_T, E, E2, alpha, D, sigma, sumK, POP)

Z=ZZ(1:end-1);
nVectI=nVect(2,:); nVectS=nVect(1,:);
nVectN=sum(nVect,1); maxN=max(nVectN);
 
I_bar=sum(nVectI.*Z')/POP;
I_T=zeros(1,length(H_T));

for i=1:length(nVectN)
    State_i = find(DemGrid==T(i));
    I_T(State_i) = I_T(State_i) + (nVectI(i)/nVectN(i))*Z(i);
end
I_T=I_T./H_T;

PQ_Sum = H_T*E;
Unstr_Inf = I_bar*alpha*sign(1-sigma)*(1-sigma)*ones(1,maxN-1); % Unstructured ext. inf.
Struc_Inf = min(sigma,1)*E*D*E2*I_T'; % Age-structured ext. inf.

% Q_ext has extra "new house" dimension
Q_ext=sparse(1,1,0,length(nVectN)+1,length(nVectN)+1);
for n=1:length(nVectN)
    if (nVect(2,n)+nVect(3,n))==0  %% as assuming infection is rare only to all S's (still in branching stage)
        R=Unstr_Inf(nVectN(n)-1)+Struc_Inf(find(DemGrid==nVectN(n)*(sumK)+nTicker(n)))/nVectN(n);
        R=R*nVectS(n); %%% This assumes external is per person %%%
        if R>0
            Q_ext(n,end)=Q_ext(n,end)+R; % Infection takes us to "new house" state
        end
    end
end

dX=(Q_int+Q_ext)'*ZZ;

end
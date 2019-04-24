function [nVect,nVectS,nVectI,nVectR,nVectN,nTicker,Q_int,Q_demo]=Add_Demography(Vect,VectS,VectI,VectR,VectN,Q_int_short, Vacc_at_birth, Birth_Ticker_Rates, Prob_at_Birth, Leave_Ticker_Rates, Prob_at_Leaving, Leave_Ticker_Rates2, Prob_at_Leaving2, Reset_Ticker_Rates, Inf_at_Reset, Ticker_at_Reset)
% Add_Demography converts the fixed-household-size internal transmission
% matrix M into a matrix nM which contains all of the transmission events
% possible for each demographic state. The associated transmission state
% vectors are converted appropriately to take total number of household
% members into account. The matrix N contains all of the demographic event
% transition rates. nTicker contains the value of the demographic ticker k
% at each state.
%   This function is called by HH_demo_structured and shouldn't need to be
% used anywhere else.

kB=size(Birth_Ticker_Rates,2); 
kL=size(Leave_Ticker_Rates,2);
kR=size(Reset_Ticker_Rates,2);

P_R=Inf_at_Reset;

Tickers=(kB)+(kL)+(kB)+(kR); % Total number of possible tick values
L=size(Vect,2); SIZE=Tickers*L; % L is number of infectious states in max size household

% First Increase the Basic Sizes
% L is size of standard system, j is Ticker position.
% So nVectX([1:L]) is the original system, nVectX(L+[1:L]) is the 2nd
% ticker etc etc.
for i=Tickers:-1:1 % Puts in full range of ticker classes for each epidemic state
    j=i-1;
    nVect(:,j*L+[1:L])=Vect;
    nVectS(j*L+[1:L])=VectS;
    nVectI(j*L+[1:L])=VectI;
    nVectR(j*L+[1:L])=VectR;
    nVectN(j*L+[1:L])=VectN;
    nTicker(j*L+[1:L])=i;
    Q_int(j*L+[1:L],j*L+[1:L])=Q_int_short; % Add a given ticker state to transmission rate matrix
end
Q_demo=0*Q_int; % Define demographic transition matrix
maxN=max(VectN);

% Vacc_at_birth -- probability vaccinated at birth is you're the first,
% second, third kid etc etc.
if length(Vacc_at_birth)==1
    Vacc_at_birth=ones(1,max(VectN))*Vacc_at_birth;
end

% Now do Births & Ticker

for i=1:(kB-1) % Again, this is for each k class to up-date the ticker
    j=i-1;
    %first remove the rate of advancement from the diagonal.
    From=[1:L]+j*L; To=[1:L]+j*L; Rate=-Birth_Ticker_Rates(VectN(1:L),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=[1:L]+j*L; To=[1:L]+j*L+L; Rate=Birth_Ticker_Rates(VectN(1:L),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end


% NOW DO IT FOR BIRTHS
i=kB; j=i-1; 
% remove the diagonals.
From=[1:L]+j*L; To=[1:L]+j*L; Rate=-Birth_Ticker_Rates(VectN(1:L),i)';
pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
From=[1:L]; 
for k=1:length(From)
    tmpVS=Vect(:,From(k)); tmpVR=tmpVS; tmpVS(1)=tmpVS(1)+1; tmpVR(3)=tmpVR(3)+1;
    Rate=Birth_Ticker_Rates(VectN(From(k)),i);
    Stop_Prob=Prob_at_Birth(end,VectN(From(k)));
    
    if VectN(From(k))==2 % stopping before any kids
        To=From(k); t=(kB)+(kL)+(kB);
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*Stop_Prob;
    else % stopping after this kid
        To=From(k); t=(kB);
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*Stop_Prob;
    end
    
    % Adding a child
    if Stop_Prob<1
        m=find(sum(abs(Vect-tmpVS*ones(1,L)),1)<0.01); % adding a susceptible
        To=m; t=0; %reset ticker
        if VectN(To)==maxN
            t=kB;
        end
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*(1-Stop_Prob)*(1-Vacc_at_birth(VectN(From(k))));
        m=find(sum(abs(Vect-tmpVR*ones(1,L)),1)<0.01); % adding a recovered/vaccinated
        To=m;
        Q_demo(From(k)+j*L,To+t*L)=Q_demo(From(k)+j*L,To+t*L)+Rate*(1-Stop_Prob)*(Vacc_at_birth(VectN(From(k))));
    end
end
        
        
% Now do First Leaving & Ticker

% FIRST CHILD LEAVING
F=find(VectN>2);
for i=1:(kL-1)
    j=i+kB-1;
    From=F+j*L; To=F+j*L; Rate=-Leave_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=F+j*L; To=F+j*L+L; Rate=Leave_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end

i=kL; j=i+kB-1;
%Now do what happens at leaving
From=find(VectN>2); clear To;
for k=1:length(From)
      Rate=Leave_Ticker_Rates(VectN(From(k)),i);
      tmpV=Vect(:,From(k));
      tmpN = sum(tmpV);
      
      if VectS(From(k))>0
          if tmpV(3)==0
              RateR=0;
          else
              if tmpV(3)==1
                  RateR=(1-P_R)/(tmpN-P_R);
              else
                  RateR=(tmpV(3)-2*P_R)/(tmpN-2*P_R);
              end
          end
          RateI = (tmpV(2)/tmpN);
          RateS = (1-RateI-RateR)*Rate; % Rate of losing a susceptible
          tmpVS = tmpV; tmpVS(1) = tmpVS(1)-1;
          m=find(sum(abs(Vect-tmpVS*ones(1,L)),1)<0.01);
          ToS = m(1);
          if VectN(From(k))==3
              t=(kB)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS,SIZE,SIZE);
          else
              t=(kB)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS,SIZE,SIZE);
          end
      end
      if VectI(From(k))>0
          RateI = (tmpV(2)/tmpN)*Rate; % Rate of losing an infected
          tmpVI = tmpV; tmpVI(2) = tmpVI(2)-1;
          m=find(sum(abs(Vect-tmpVI*ones(1,L)),1)<0.01);
          ToI = m(1);
          if VectN(From(k))==3
              t=(kB)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI,SIZE,SIZE);
          else
              t=(kB)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI,SIZE,SIZE);
          end
      end
      if VectR(From(k))>0
          if tmpV(3)==0
              RateR=0;
          else
              if tmpV(3)==1
                  RateR=Rate*(1-P_R)/(tmpN-P_R);
              else
                  RateR=Rate*(tmpV(3)-2*P_R)/(tmpN-2*P_R);
              end
          end
          tmpVR = tmpV; tmpVR(3) = tmpVR(3)-1;
          m=find(sum(abs(Vect-tmpVR*ones(1,L)),1)<0.01);
          ToR = m(1);
          if VectN(From(k))==3
              t=(kB)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR,SIZE,SIZE);
          else
              t=(kB)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR,SIZE,SIZE);
          end
      end
      
end
      
 % NOW OTHER CHILDREN LEAVING     
               
F=find(VectN>2);
for i=1:(kB-1)
    j=i+kB+kL-1;
    From=F+j*L; To=F+j*L; Rate=-Leave_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=F+j*L; To=F+j*L+L; Rate=Leave_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end

i=kB; j=i+kB+kL-1;
%Now do what happens at leaving
From=find(VectN>2); clear To;
for k=1:length(From)
      Rate=Leave_Ticker_Rates(VectN(From(k)),i);
      tmpV=Vect(:,From(k));
      tmpN = sum(tmpV);
      
      if VectS(From(k))>0
          if tmpV(3)==0
              RateR=0;
          else
              if tmpV(3)==1
                  RateR=(1-P_R)/(tmpN-P_R);
              else
                  RateR=(tmpV(3)-2*P_R)/(tmpN-2*P_R);
              end
          end
          RateI = (tmpV(2)/tmpN);
          RateS = (1-RateI-RateR)*Rate; % Rate of losing a susceptible
          tmpVS = tmpV; tmpVS(1) = tmpVS(1)-1;
          m=find(sum(abs(Vect-tmpVS*ones(1,L)),1)<0.01);
          ToS = m(1);
          if VectN(From(k))==3
              t=(kB)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS,SIZE,SIZE);
          else
              t=(kB)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToS+t*L,RateS,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateS,SIZE,SIZE);
          end
      end
      if VectI(From(k))>0
          RateI = (tmpV(2)/tmpN)*Rate; % Rate of losing an infected
          tmpVI = tmpV; tmpVI(2) = tmpVI(2)-1;
          m=find(sum(abs(Vect-tmpVI*ones(1,L)),1)<0.01);
          ToI = m(1);
          if VectN(From(k))==3
              t=(kB)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI,SIZE,SIZE);
          else
              t=(kB)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToI+t*L,RateI,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateI,SIZE,SIZE);
          end
      end
      if VectR(From(k))>0
          if tmpV(3)==0
              RateR=0;
          else
              if tmpV(3)==1
                  RateR=Rate*(1-P_R)/(tmpN-P_R);
              else
                  RateR=Rate*(tmpV(3)-2*P_R)/(tmpN-2*P_R);
              end
          end
          tmpVR = tmpV; tmpVR(3) = tmpVR(3)-1;
          m=find(sum(abs(Vect-tmpVR*ones(1,L)),1)<0.01);
          ToR = m(1);
          if VectN(From(k))==3
              t=(kB)+(kL)+(kB);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR,SIZE,SIZE);
          else
              t=(kB)+(kL);
              Q_demo=Q_demo+sparse(From(k)+j*L,ToR+t*L,RateR,SIZE,SIZE); Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-RateR,SIZE,SIZE);
          end
      end
      
end

       
% NOW DO RESET AND TICKER
F=find(VectN==2); % must be just 2 adults.
for i=1:(kR-1)
    j=i+kB+kL+kB-1;
    From=F+j*L; To=F+j*L; Rate=-Reset_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
    
    From=F+j*L; To=F+j*L+L; Rate=Reset_Ticker_Rates(VectN(F),i)';
    pos=find(From>0 & To>0 & Rate~=0); Q_demo=Q_demo+sparse(From(pos),To(pos),Rate(pos),SIZE,SIZE);
end

%Now do what happens at reset
i=kR; j=i+kB+kL+kB-1;
From=find(VectN==2);
ToProb(3)=(1-Inf_at_Reset).^2; To(3)=0; % both sus
m=find(VectS==2 & VectI==0 & VectR==0); % Passed over by outbreak
if ~isempty(m)
    To(3)=m;
end

ToProb(2)=2*(Inf_at_Reset)*(1-Inf_at_Reset);
m=find(VectS==1 & VectI==0 & VectR==1); % 1 recovered
if ~isempty(m)
    To(2)=m;
end

ToProb(1)=(Inf_at_Reset).^2;
m=find(VectS==0 & VectI==0 & VectR==2); % 2 recovers
if ~isempty(m)
    To(1)=m;
end

Rate=Reset_Ticker_Rates(VectN(From),i);
for k=1:length(From)
    Q_demo=Q_demo+sparse(From(k)+j*L,From(k)+j*L,-Rate(k),SIZE,SIZE);
    for K=1:3
       if Rate(k)*ToProb(K)>0 
           Q_demo=Q_demo+sparse(From(k)+j*L,To(K),Rate(k)*ToProb(K),SIZE,SIZE);% N=N+sparse(From(k),From(k),-Rate(k)*ToProb(K),SIZE,SIZE);
       end
    end
end
       
% NOW NEED TO REMOVE SURPLUS POINTS
% if waiting for kids to leave and only size 2
r=find(nVectN==2 & (nTicker>kB & nTicker<=kB+kL+kB));
% if waiting for second kids to leave and size=max
r=[r find(nVectN==max(nVectN) & (nTicker>kB+kL & nTicker<=kB+kL+kB))];
% if waiting for birth and size = max
r=[r find(nVectN==max(nVectN) & nTicker<=kB)];
% if waiting for reset and bigger than size 2
r=[r find(nVectN>2 & nTicker>kB+kL+kB)];

nVect(:,r)=[]; % This syntax just deletes these entries (changing size of array)  
nVectS(r)=[]; 
nVectI(r)=[]; 
nVectR(r)=[]; 
nVectN(r)=[]; 
nTicker(r)=[];
Q_int(r,:)=[]; Q_int(:,r)=[];
Q_demo(r,:)=[]; Q_demo(:,r)=[];

end

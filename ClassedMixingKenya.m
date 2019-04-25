% This script calculates the demographic class-stratified transmission
% matrix E, along with the proportion of individuals in each age class, for
% the Kenya-like demographic setting.

ModelBds=[0 1 6 16 20 50 99]; % Age class boundaries from Kilifi paper

% We begin by calculating the age-structured WAIFW matrix
load('ContactData/KenyaContactData.mat');
good=find(MyID~=MyID(find(isnan(SameHH)==1)))';
ContactAge=ContactAge(good);
MyAge=MyAge(good)';
MyID=MyID(good)';
SameHH=SameHH(good);
Ext_Contact=1-SameHH;
D_Ext=zeros(length(ModelBds)-1,length(ModelBds)-1); % WAIFW matrix with only external contacts
D_Int=D_Ext; % Only internal contacts (needed for d_int)
D_All=D_Ext; % This includes all interactions
Contact_No_Ext=D_Ext;
Contact_No_Int=D_Int;
Contact_No_All=D_Ext;

Contact_Duration=[];
Contact_Duration(Ext_Contact==0)=0.1018;
Contact_Duration(Ext_Contact==1)=0.0577;
Total_Duration_internal=full(sparse(MyID(Ext_Contact==0),1,Contact_Duration(Ext_Contact==0),max(MyID),1));
Total_Duration_external=full(sparse(MyID(Ext_Contact==1),1,Contact_Duration(Ext_Contact==1),max(MyID),1));

d_int=mean(Total_Duration_internal(unique(MyID)));
d_ext=mean(Total_Duration_external(unique(MyID)));
d_all=mean(Contact_Duration); % We need this to do purely age-structured mixing

disp('Now calculating int strength matrix')
for Class1=1:length(ModelBds)-1
    Participant_in_C1=MyAge==Class1; % Find contact events where participant is in class
    for Class2=1:length(ModelBds)-1
        Contact_in_C2=ContactAge==Class2; % Find contact events where contact is in class
        D_Ext(Class1,Class2)=(Participant_in_C1.*Contact_in_C2.*Ext_Contact)*Contact_Duration';
        D_Int(Class1,Class2)=(Participant_in_C1.*Contact_in_C2.*(1-Ext_Contact))*Contact_Duration';
        D_All(Class1,Class2)=(Participant_in_C1.*Contact_in_C2)*Contact_Duration';
        
    end
    Total_in_C1=length(unique(MyID(Participant_in_C1))); % This is total number of participants in class
    D_Ext(Class1,:)=D_Ext(Class1,:)/Total_in_C1; % Scale by number of contacts recorded for this age class
    D_Int(Class1,:)=D_Int(Class1,:)/Total_in_C1;
    D_All(Class1,:)=D_All(Class1,:)/Total_in_C1;
end

disp('int strength matrix calculated')

TB=(35.2/12)*365;
kB=2;
TL=22.7*365;
kL=4;
TD=72.45*365;
kR=5;
Distrib_Children=[0.019 0.038 0.098 0.129 0.15 0.132 0.119 0.108 0.078 0.05 0.079];
Constrained_max=floor(TL/TB)+1; % This is the most kids we can have under our model assumptions with these parameter choices
Distrib_Children(Constrained_max+1)=sum(Distrib_Children(Constrained_max+1:end)); % +1 since first one is no kids
Distrib_Children=Distrib_Children(1:Constrained_max+1);
Distrib_Children=Distrib_Children/sum(Distrib_Children);
Cond_DC=Distrib_Children(2:end)/sum(Distrib_Children(2:end)); % Distribution of number of children conditioned on having kids
tmp=[0 cumsum(Distrib_Children)];
StopProb=[0 Distrib_Children./(1-tmp(1:(end-1)))]; % Prob of stopping after each number of kids
StopProb(end)=1;

disp('Now calculating translation matrix')
points=100000; % This is number of points for Monte Carlo integration
[E,NGrid,tickGrid,DemGrid] = State_to_Class_Mat_With_Elders(kB,kL,kR,TB,TL,365*ModelBds,StopProb,points);
disp('Translation matrix calculated')

disp('Now calculating age class proportions')
ClassProb=zeros(1,length(ModelBds)-1); % This loop will do age class proportions for children
for i=1:length(ModelBds)-2 % minus 2 since we assume no kids in elder class
    tic
    L=365*ModelBds(i); U=365*ModelBds(i+1);
    fun1 = @(A,L,kB,kL,TB,TL,Cond_DC) LeaveAgeProb(A,kB,kL,TB,TL,Cond_DC).*((A-L)./A); % Integrand for leaving within age class
    fun2 = @(A,L,U,kB,kL,TB,TL,Cond_DC) LeaveAgeProb(A,kB,kL,TB,TL,Cond_DC).*((U-L)./A); % Integrand for leaving after age class
    int1=integral(@(A)fun1(A,L,kB,kL,TB,TL,Cond_DC),L,U);
    int2=integral(@(A)fun2(A,L,U,kB,kL,TB,TL,Cond_DC),U,TD); % Assume people have definitely left home by LE
    ClassProb(i)=int1+int2;
    toc
end
% Now assume adults spend TR-TL alive after leaving home
ClassProb=ClassProb*TL/TD; % Scale by proportion of lifetime spent as child
ClassProb(end-1)=ClassProb(end-1)+(365*ModelBds(end-1)-TL)/TD;
ClassProb(end)= 1-sum(ClassProb);
disp('Age class proportions calculated')

filename=['Parameters/Kenya_MixingData_' datestr(now,'ddmmyy_HHMMSS')];
save(filename, 'D_Ext', 'D_All', 'd_ext',...
    'd_int', 'd_all', 'ClassProb', 'E', 'NGrid', 'tickGrid', 'DemGrid');

function [P_Age] = LeaveAgeProb(A,kB,kL,TB,TL,Cond_DC)
% Calculates probability of kid leaving home at age A by going over all
% possible numbers of children and positions
% Recall that Cond_DC is conditioned on having at least one child!
% Nmax is max number of children and N is number of children. n is position
% of current child in birth ordering.
Nmax=length(Cond_DC);
P_Age=0;
for N=1:Nmax
    summand=0;
    for n=1:N
        lam=[(kB/TB)*ones(1,(N-n)*kB) (kL/(TL-(N-1)*TB))*ones(1,kL) (kB/TB)*ones(1,(n-1))*kB];
        summand=summand+(1/N)*HypoExpPdf(A,lam);
    end
    P_Age=P_Age+Cond_DC(N)*summand;
end
end

function [P_Hyp] = HypoExpPdf(x,lambda)
% Hypoexponential distribution - see wikipedia etc
range = length(x);
n = length(lambda);

alpha = zeros(1,n);
alpha(1) = 1;

index = 1:n;
Theta = sparse(index,index,-lambda,n,n);
Theta = Theta + sparse(index(1:end-1),index(1:end-1)+1,lambda(1:end-1),n,n);

for l=range:-1:1
    xTheta = x(l)*Theta;
    xTheta = expm(xTheta);
    P_Hyp(l) = -alpha*xTheta*Theta*ones(n,1);
end

P_Hyp = full(P_Hyp);

end
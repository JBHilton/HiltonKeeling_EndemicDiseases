% This script calculates the demographic class-stratified transmission
% matrix E, along with the proportion of individuals in each age class, for
% the UK-like demographic setting.

ModelBds=[0 1 2 4 10 17 98]; % Age class boundaries for children

TRANSLATE_CONTACT=[2.5 10 30 120 240]; % This translates the contact duration categories from POLYMOD to a duration in minutes

% We begin by calculating the age-structured WAIFW matrix
load('ContactData/Polymod_data.mat');
good=find(Contact_Age==Contact_Age&Contact_Duration==Contact_Duration...
    &Contact_Freq==Contact_Freq&Contact_Loc_Main==Contact_Loc_Main&My_Age==My_Age); % This gets rid of any missing data
Global_ID=Global_ID(good);
Contact_Age=Contact_Age(good);
Contact_Duration=Contact_Duration(good);
Contact_Freq=Contact_Freq(good);
My_Age=My_Age(good);
Contact_Loc_Main=Contact_Loc_Main(good);
Ext_Contact=(Contact_Loc_Main~=1); % Finds all the contacts which are not home contacts
D_Ext_Unscaled=zeros(length(ModelBds)-1,length(ModelBds)-1); % WAIFW matrix with only external contacts
D_Int_Unscaled=D_Ext_Unscaled; % Only internal contacts (needed to get duration for internal mixing)
D_All_Unscaled=D_Ext_Unscaled; % This includes all interactions
Contact_No_Ext=D_Ext_Unscaled;
Contact_No_Int=D_Int_Unscaled;
Contact_No_All=D_Ext_Unscaled;

Contact_Duration=TRANSLATE_CONTACT(Contact_Duration)/(24*60); % Rescale to units of days

% Note that sparse sums over repeat indices so this does sum of all a
% participants internal and external contacts
Total_Duration_internal=full(sparse(Global_ID(Ext_Contact==0),1,Contact_Duration(Ext_Contact==0),max(Global_ID),1));
Total_Duration_external=full(sparse(Global_ID(Ext_Contact==1),1,Contact_Duration(Ext_Contact==1),max(Global_ID),1));
Total_Duration_all=full(sparse(Global_ID,1,Contact_Duration,max(Global_ID),1));

d_int=mean(Total_Duration_internal);
d_ext=mean(Total_Duration_external);
d_all=mean(Total_Duration_all); % We need this to do purely age-structured mixing
d_int_kids=mean(Contact_Duration(Ext_Contact==0&My_Age<ModelBds(end-1)));% Needed for tau calculation, this conditions on being in lower age classes

% This loop calculates WAIFW matrices
disp('Now calculating duration matrix')
for Class1=1:length(ModelBds)-1 % Indexed by age class
    Participant_in_C1=My_Age>=ModelBds(Class1)&My_Age<ModelBds(Class1+1); % Find contact events where participant is in class
    for Class2=1:length(ModelBds)-1
        Contact_in_C2=Contact_Age>=ModelBds(Class2)&Contact_Age<ModelBds(Class2+1); % Find contact events where contact is in class
        D_Ext_Unscaled(Class1,Class2)=(Participant_in_C1.*Contact_in_C2.*Ext_Contact)*Contact_Duration';
        D_Int_Unscaled(Class1,Class2)=(Participant_in_C1.*Contact_in_C2.*(1-Ext_Contact))*Contact_Duration';
        D_All_Unscaled(Class1,Class2)=(Participant_in_C1.*Contact_in_C2)*Contact_Duration';
        
    end
    Total_in_C1=length(unique(Global_ID(Participant_in_C1))); % This is total number of participants in class
    D_Ext_Unscaled(Class1,:)=D_Ext_Unscaled(Class1,:)/Total_in_C1; % Scale by number of contacts recorded for this age class
    D_Int_Unscaled(Class1,:)=D_Int_Unscaled(Class1,:)/Total_in_C1;
    D_All_Unscaled(Class1,:)=D_All_Unscaled(Class1,:)/Total_in_C1;
end

disp('Duration matrix calculated')

% UK-like demographic parameters
TB=3*365;
kB=2;
TL=25*365;
kL=4;
TD=80*365;
kR=5;
Distrib_Children=[0.17 0.18 0.37 0.17 0.1];
Distrib_Children=Distrib_Children/sum(Distrib_Children);
Cond_DC=Distrib_Children(2:end)/sum(Distrib_Children(2:end)); % Distribution of number of children conditioned on having kids
tmp=[0 cumsum(Distrib_Children)];
StopProb=[0 Distrib_Children./(1-tmp(1:(end-1)))]; % Prob of stopping after each number of kids
StopProb(end)=1;

disp('Now calculating translation matrix')
points=100000; % This is number of points for Monte Carlo integration
[E,NGrid,tickGrid,DemGrid] = State_to_Class_Mat(kB,kL,kR,TB,TL,365*ModelBds,StopProb,points); % This calculates the translation matrix E
disp('Translation matrix calculated')

disp('Now calculating age class proportions')
ClassProb=zeros(1,length(ModelBds)-1); % This loop will do age class proportions for children
for i=1:length(ModelBds)-1
    tic
    L=365*ModelBds(i); U=365*ModelBds(i+1);
    fun1 = @(A,L,kB,kL,TB,TL,Cond_DC) LeaveAgeProb(A,kB,kL,TB,TL,Cond_DC).*((A-L)./A); % Integrand for leaving within age class
    fun2 = @(A,L,U,kB,kL,TB,TL,Cond_DC) LeaveAgeProb(A,kB,kL,TB,TL,Cond_DC).*((U-L)./A); % Integrand for leaving after age class
    int1=integral(@(A)fun1(A,L,kB,kL,TB,TL,Cond_DC),L,U);
    int2=integral(@(A)fun2(A,L,U,kB,kL,TB,TL,Cond_DC),U,150*365); % Assume people have definitely left home by LE
    ClassProb(i)=int1+int2;
    toc
end
% Now assume adults spend TR-TL alive after leaving home
ClassProb=ClassProb*TL/TD; % Scale by proportion of lifetime spent as child
ClassProb(end)=ClassProb(end)+(TD-TL)/TD;
disp('Age class proportions calculated')

filename=['UK_MixingData_' datestr(now,'ddmmyy_HHMMSS')];
save(filename, 'D_Ext_Unscaled', 'D_All_Unscaled', 'd_ext',...
    'd_int', 'd_int_kids', 'd_all', 'ClassProb', 'E',...
    'NGrid', 'tickGrid', 'DemGrid');

function [P_Age] = LeaveAgeProb(A,kB,kL,TB,TL,Cond_DC)
% Calculates probability of kid leaving home at age A by going over all
% possible numbers of children and positions
% Recall that Cond_DC is conditioned on having at least one child!
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
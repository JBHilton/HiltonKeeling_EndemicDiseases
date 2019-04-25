function [ b_young,b_old,fig ] = Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil,DiseaseFree,figname,pltname,varargin )
%Get_Prevalence_Histogram creates a stacked histogram of infectious
%prevalence given the equilibrium distribution Equil and the disease-free
%equilibrium distribution DiseaseFree.
%   Set LEGEND=1 to plot legend on histogram, or LEGEND=0 to leave it off.
%   figname is a character string which tells us what to save figure etc
%   as. Set CONDITIONAL=1 to plot distributions conditioned on household
%   size distribution, CONDITIONAL=0 to plot raw percentages

ip=inputParser;

defaultYoungMax=[];
defaultOldMax=[];

defaultLegend='off';
validLegend={'on','off'};
checkLegend=@(x) any(validatestring(x,validLegend));

defaultCond='cond';
validCond={'cond','abs'};
checkCond=@(x) any(validatestring(x,validCond));

defaultYmax=1;

defaultYOutline=[];
defaultOOutline=[];

defaultSave='Off';
validSave={'On','Off'};
checkSave=@(x) any(validatestring(x,validSave));

addRequired(ip,'nVect',@isnumeric);
addRequired(ip,'nTicker',@isnumeric);
addRequired(ip,'kB',@isnumeric);
addRequired(ip,'kL',@isnumeric);
addRequired(ip,'Equil',@isnumeric);
addRequired(ip,'DiseaseFree',@isnumeric);
addRequired(ip,'figname',@ischar);
addRequired(ip,'pltname',@ischar);
addOptional(ip,'young_max',defaultYoungMax,@isnumeric);
addOptional(ip,'old_max',defaultOldMax,@isnumeric);
addOptional(ip,'LEGEND',defaultLegend,checkLegend);
addOptional(ip,'CONDITIONAL',defaultCond,checkCond);
addOptional(ip,'y_max',defaultYmax,@isnumeric);
addOptional(ip,'YOutline',defaultYOutline);
addOptional(ip,'OOutline',defaultOOutline);
addOptional(ip,'Save',defaultSave,checkSave);

parse(ip,nVect,nTicker,kB,kL,Equil,DiseaseFree,figname,pltname,varargin{:});
young_max=ip.Results.young_max; old_max=ip.Results.old_max;
LEGEND=ip.Results.LEGEND; CONDITIONAL=ip.Results.CONDITIONAL;
y_max=ip.Results.y_max;
YOutline=ip.Results.YOutline; OOutline=ip.Results.OOutline;
Save=ip.Results.Save;

ColourMap1=jet; ColourMap1=ColourMap1(1:7:end,:);
ColourMap2=jet; ColourMap2=ColourMap2(end:-7:1,:);

nVectN=sum(nVect,1);
nVectI=nVect(2,:);

maxN=max(nVectN);

Index=cell(1,2*maxN-3);
x_lab=cell(1,2*maxN-3);
Index{1}=find(nVectN==2&nTicker<kB);
x_lab{1}=num2str(2);
Index{2*maxN-3}=find(nVectN==2&nTicker>=kB);
x_lab{2*maxN-3}=num2str(2);
Index{maxN-1}=find(nVectN==maxN);
x_lab{maxN-1}=num2str(maxN);
for i=3:maxN-1
    Index{i-1}=find(nVectN==i&nTicker<kB+kL);
    x_lab{i-1}=num2str(i);
    Index{2*(maxN)-i-1}=find(nVectN==i&nTicker>=kB+kL);
    x_lab{2*(maxN)-i-1}=num2str(i);
end
if strcmp(CONDITIONAL,'cond')
    y_DiseaseFree=zeros(1,2*maxN-3);
    for i=1:2*maxN-3
        y_DiseaseFree(i)=sum(DiseaseFree(Index{i}));
    end
end

y_Equil_young=zeros(maxN-1,maxN);
y_Equil_old=zeros(maxN-2,maxN);
Abs_Freq=zeros(1,2*maxN-3);
for i=1:maxN-1
    for j=1:i+1
        if strcmp(CONDITIONAL,'cond')
            y_Equil_young(i,j)=100*sum(Equil(Index{i}(nVectI(Index{i})==j)))/y_DiseaseFree(i);
        else
            y_Equil_young(i,j)=100*sum(Equil(Index{i}(nVectI(Index{i})==j)));
        end
        Abs_Freq(i)=Abs_Freq(i)+y_Equil_young(i,j)*j;
    end
end
for i=1:maxN-2
    for j=1:i+1
        if strcmp(CONDITIONAL,'cond')
            y_Equil_old(i,j)=100*sum(Equil(Index{maxN-1+i}(nVectI(Index{maxN-1+i})==j)))/y_DiseaseFree(maxN-1+i);
        else
            y_Equil_old(i,j)=100*sum(Equil(Index{maxN-1+i}(nVectI(Index{maxN-1+i})==j)));
        end
        Abs_Freq(maxN-1+i)=Abs_Freq(maxN-1+i)+y_Equil_old(i,j)*j;
    end
end

if isempty(young_max)
    young_max=find(max(y_Equil_young)<1e-3,1)-1;
end
if isempty(old_max)
    old_max=find(max(y_Equil_old)<1e-3,1)-1;
end

y_Equil_young=y_Equil_young(:,1:young_max); y_Equil_old=y_Equil_old(:,1:old_max);
 
fig=figure('DefaultAxesFontSize',24,'units','normalized','outerposition',[0 0 1 1]);
b_young=bar(1:maxN-1,y_Equil_young,'stacked');
hold on
b_old=bar(maxN:2*maxN-3,y_Equil_old,'stacked');
bar_label=strings(1,young_max+old_max);
for i=1:young_max % This for loop does the colours
    set(b_young(i),'FaceColor',ColourMap1(i,:));
    bar_label(i)=[num2str(i) ' Infected'];
end
for i=1:old_max
    set(b_old(i),'FaceColor',ColourMap2(i,:));
    bar_label(young_max+i)=[num2str(i) ' Infected'];
end
line([maxN-1,maxN-1],[0,1],'Color','k','LineStyle','--','LineWidth',1) % Divider between young and old households
plot(Abs_Freq,'ko','LineWidth',2,'MarkerFaceColor','w','MarkerSize',10);
if strcmp(LEGEND,'on')
    l=legend([b_young b_old],bar_label,'Location','northeast');
end
if ~isempty(OOutline)
    ydat=zeros(1,maxN-2);
    for stack=1:size(OOutline,2)
        ydat=ydat+stack*OOutline(stack).YData;
    end
    bar(OOutline(1).XData,ydat,1,'LineWidth',2,'EdgeColor','m','FaceColor','none');
end
if ~isempty(YOutline)
    ydat=zeros(1,maxN-1);
    for stack=1:size(YOutline,2)
        ydat=ydat+stack*YOutline(stack).YData;
    end
    bar(YOutline(1).XData,ydat,1,'LineWidth',2,'EdgeColor','m','FaceColor','none');
end
axis([0 2*maxN-2 0 y_max])
axis square
xticks(1:2*maxN-3)
xticklabels(x_lab)
xlabel('Household size')
if strcmp(CONDITIONAL,'cond')
    ylabel('Conditional Frequency (%)');
else
    ylabel('Absolute Frequency (%)');
end
ylim=get(gca,'ylim');
if strcmp(LEGEND,'off')
    t=text(maxN-1,ylim(2),pltname,'FontSize',24);
    pause(1); % This is here so figure loads properly before we do anything complicated with positions
    t.Position(1:2)=t.Position(1:2)-0.5*t.Extent(3:4);
else
    t=text(maxN-1,ylim(2),pltname,'FontSize',24);
    pause(1); % This is here so figure loads properly before we do anything complicated with positions
    t.Position(1)=t.Position(1)-1.1*t.Extent(3);
    t.Position(2)=t.Position(2)-0.5*t.Extent(4);
end

if strcmp(Save,'On')
    savefig(fig, figname);
    fig.PaperPositionMode = 'auto'; % This bit gets rid of whitespace
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,figname,'-depsc');
end

end


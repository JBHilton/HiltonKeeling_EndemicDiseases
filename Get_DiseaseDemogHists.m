% This script plots the equilibrium disease burden histograms for measles
% and mumps in the UK and Kenya

% The homogeneous histograms are calculated to provide the 'pink bars'
% comparisons in the paper

% First do UK
load('Parameters/UK_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/UKDemography.mat')

load('ModelOutput/UKMeaslesEqDist.mat')
load('ModelOutput/UKMeaslesHomEqDist.mat');
[UKMeasles_youngC,UKMeasles_oldC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'UKMeaslesHist_HomC','Homogeneous Mixing','y_max',0.8 );
[~,~,h_UKMeasles_C]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MeaslesUK,DiseaseFree,'UKMeaslesHist_Cond','UK, Measles','y_max',0.8,'YOutline',UKMeasles_youngC,'OOutline',UKMeasles_oldC );
[UKMeasles_youngUC,UKMeasles_oldUC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'UKMeaslesHist_HomUC','Homogeneous Mixing','y_max',0.03,'CONDITIONAL','abs' );
[~,~,h_UKMeasles_UC]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MeaslesUK,DiseaseFree,'UKMeaslesHist_Unco','UK, Measles','y_max',0.03,'CONDITIONAL','abs','YOutline',UKMeasles_youngUC,'OOutline',UKMeasles_oldUC );

load('ModelOutput/UKMumpsEqDist.mat')
load('ModelOutput/UKStructureEqDists.mat');
[UKMumps_youngC,UKMumps_oldC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'UKMumpsHist_HomC','Homogeneous Mixing','y_max',0.8 );
[~,~,h_UKMumps_C]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MumpsUK,DiseaseFree,'UKMumpsHist_Cond','UK, Mumps','y_max',0.8,'YOutline',UKMumps_youngC,'OOutline',UKMumps_oldC );
[UKMumps_youngUC,UKMumps_oldUC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'UKMumpsHist_HomUC','Homogeneous Mixing','y_max',0.03,'CONDITIONAL','abs' );
[~,~,h_UKMumps_UC]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MumpsUK,DiseaseFree,'UKMumpsHist_Unco','UK, Mumps','y_max',0.03,'CONDITIONAL','abs','YOutline',UKMumps_youngUC,'OOutline',UKMumps_oldUC );

%%
% Now Kenya

load('Parameters/Kenya_MixingData.mat'); % Contains d_ext, d_int, ClassProb, D_All, D_Ext, E, NGrid, tickGrid, DemGrid

load('Parameters/KenyaDemography.mat')

load('ModelOutput/KenyaMeaslesEqDist.mat')
load('ModelOutput/KenyaMeaslesHomEqDist.mat');
[KenyaMeasles_youngC,KenyaMeasles_oldC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'KenyaMeaslesHist_HomC','Homogeneous Mixing','y_max',0.8 );
[~,~,h_KenyaMeasles_C]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MeaslesKenya,DiseaseFree,'KenyaMeaslesHist_Cond','Kenya, Measles','y_max',0.8,'LEGEND','on','young_max',4,'old_max',3,'YOutline',KenyaMeasles_youngC,'OOutline',KenyaMeasles_oldC );
[KenyaMeasles_youngUC,KenyaMeasles_oldUC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'KenyaMeaslesHist_HomUC','Homogeneous Mixing','y_max',0.03,'CONDITIONAL','abs' );
[~,~,h_KenyaMeasles_UC]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MeaslesKenya,DiseaseFree,'KenyaMeaslesHist_Unco','Kenya, Measles','y_max',0.03,'LEGEND','on','young_max',4,'old_max',2,'CONDITIONAL','abs','YOutline',KenyaMeasles_youngUC,'OOutline',KenyaMeasles_oldUC );

load('ModelOutput/KenyaMumpsEqDist.mat')
load('ModelOutput/KenyaMumpsHomEqDist.mat');
[KenyaMumps_youngC,KenyaMumps_oldC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'KenyaMumpsHist_HomC','Homogeneous Mixing','y_max',0.8 );
[~,~,h_KenyaMumps_C]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MumpsKenya,DiseaseFree,'KenyaMumpsHist_Cond','Kenya, Mumps','y_max',0.8,'YOutline',KenyaMumps_youngC,'OOutline',KenyaMumps_oldC );
[KenyaMumps_youngUC,KenyaMumps_oldUC,~]=Get_Prevalence_Histogram...
    ( nVect,nTicker,kB,kL,Equil_Hom,DiseaseFree,'KenyaMumpsHist_HomUC','Homogeneous Mixing','y_max',0.03,'CONDITIONAL','abs' );
[~,~,h_KenyaMumps_UC]=Get_Prevalence_Histogram( nVect,nTicker,kB,kL,Equil_MumpsKenya,DiseaseFree,'KenyaMumpsHist_Unco','Kenya, Mumps','y_max',0.03,'CONDITIONAL','abs','YOutline',KenyaMumps_youngUC,'OOutline',KenyaMumps_oldUC );
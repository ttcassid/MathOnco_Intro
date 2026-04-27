%% This script plots the tumour size dynamics from Craig et al. PLOS Comp Biol 2019
close all
clear all
tic
% Colours for plotting
ColorMatrix = [141,211,199;
    253,180,98;
    190,186,218 ;
    251,128,114 ;
    128,177,211;
    178,223,138 ;
    51,160,44 ;
    227,26,28;
    255,127,0 ;
    106,61,154]./255;
%% Load data for monoculture data
MutantSwitch = 1; % 0 for WT, 1 for M1, 2 for M2, 3 for M3
TreatmentColorIndex = 1;
TimeVec = [1,3,5,7];
NDataPoints = length(TimeVec);

% Untreated data
UntreatedDataMatrix = struct2array(load('CAT_Monoculture_NoTreatment_Data.mat'));
UntreatedDataMutant = UntreatedDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(UntreatedDataMatrix(2+MutantSwitch,1))); % Normalized to fold change

% Docetaxel data
DocetaxelDataMatrix = struct2array(load('CAT_Monoculture_Docetaxel_Data.mat'));
DocetaxelDataMutant = DocetaxelDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(DocetaxelDataMatrix(2+MutantSwitch,1))); % Normalized to fold change

%Afatinib
AfatinibDataMatrix = struct2array(load('CAT_Monoculture_Afatinib_Data.mat'));
AfatinibDataMutant = AfatinibDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(AfatinibDataMatrix(2+MutantSwitch,1))); % Normalized to fold change
%%
Fig1 =  figure(1);
g60 = scatter(TimeVec,UntreatedDataMutant,60,ColorMatrix(TreatmentColorIndex, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
g60 = scatter(TimeVec(NDataPoints-1:end),DocetaxelDataMutant(NDataPoints-1:end),60,ColorMatrix(TreatmentColorIndex + 1, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
g60 = scatter(TimeVec(NDataPoints-1:end),AfatinibDataMutant(NDataPoints-1:end),60,ColorMatrix(TreatmentColorIndex + 2, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
ylim([10^2.5 10^5])
set(gca,'yscale','log')





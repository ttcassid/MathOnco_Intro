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
UntreatedDataMatrix = struct2array(load('CAT_50WT_50M1_Data.mat'));

UntreatedDataCoculture_WT = UntreatedDataMatrix(2,:).*(UntreatedDataMatrix(2,1)./(UntreatedDataMatrix(2,1))); % Normalized to fold change
UntreatedDataCoculture_Mutant = UntreatedDataMatrix(3,:).*(UntreatedDataMatrix(3,1)./(UntreatedDataMatrix(3,1))); % Normalized to fold change
UntreatedTotal = UntreatedDataCoculture_WT + UntreatedDataCoculture_Mutant;

DocetaxelDataCoculture_WT = UntreatedDataMatrix(4,:).*(UntreatedDataMatrix(2,1)./(UntreatedDataMatrix(4,1))); % Normalized to fold change
DocetaxelDataCoculture_Mutant = UntreatedDataMatrix(5,:).*(UntreatedDataMatrix(3,1)./(UntreatedDataMatrix(5,1))); % Normalized to fold change
DocetaxelTotal = DocetaxelDataCoculture_WT + DocetaxelDataCoculture_Mutant;

AfatinibDataCoculture_WT = UntreatedDataMatrix(6,:).*(UntreatedDataMatrix(2,1)./(UntreatedDataMatrix(6,1))); % Normalized to fold change
AfatinibDataCoculture_Mutant = UntreatedDataMatrix(7,:).*(UntreatedDataMatrix(3,1)./(UntreatedDataMatrix(7,1))); % Normalized to fold change
AfatinibTotal = AfatinibDataCoculture_WT + AfatinibDataCoculture_Mutant;

%%
Fig1 =  figure(1);
g60 = scatter(TimeVec,UntreatedTotal,60,ColorMatrix(TreatmentColorIndex, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
g1 = scatter(TimeVec,UntreatedDataCoculture_WT,40,ColorMatrix(10, :),'s','LineWidth',1.85);
hold on
g1 = scatter(TimeVec,UntreatedDataCoculture_Mutant,40,ColorMatrix(9, :),'s','LineWidth',1.85);
hold on
ylim([10^1.5 10^5])
set(gca,'yscale','log');

Fig2 =  figure(2);
g60 = scatter(TimeVec,DocetaxelTotal,60,ColorMatrix(TreatmentColorIndex + 1, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on 
g1 = scatter(TimeVec,DocetaxelDataCoculture_WT,40,ColorMatrix(10, :),'s','LineWidth',1.85);
hold on
g1 = scatter(TimeVec,DocetaxelDataCoculture_Mutant,40,ColorMatrix(9, :),'s','LineWidth',1.85);
hold on
ylim([10^1.5 10^5]);
set(gca,'yscale','log');

Fig3 = figure(3);
g60 = scatter(TimeVec,AfatinibTotal,60,ColorMatrix(TreatmentColorIndex + 2, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
g1 = scatter(TimeVec,AfatinibDataCoculture_WT,40,ColorMatrix(10, :),'s','LineWidth',1.85);
hold on
g1 = scatter(TimeVec,AfatinibDataCoculture_Mutant,40,ColorMatrix(9, :),'s','LineWidth',1.85);
hold on
ylim([10^1.5 10^5])
set(gca,'yscale','log')





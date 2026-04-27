%% This script plots the dose response data from Kavran et al. PNAS 2022
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
% MutantSwitch = 1; % 0 for WT, 1 for M1, 2 for M2, 3 for M3
TreatmentColorIndex = 1; % Choose the colour
DoseVec = [0.01,0.1,0.3,1,3,10,30,100,300,1000,3000,10000]; % Dose sizes
NDataPoints = length(DoseVec); %How many distinct doses considered

% MEKi data
DoseResponseDataMatrix = struct2array(load('BRAFi_DoseResponse_Mutant.mat')); % table2array(struct2table(load('BRAFi_DoseResponse_Mutant.mat')))
DoseResponseData = DoseResponseDataMatrix; %load the dose response data
DoseVec = DoseResponseData(1,2:end);
NReplicates = 4; % how many replicates are there 
for ii = 1:1; %  6; % 6 total cell lines
    CellLine = ii; % 2 ;  % 1 : WT, 2 :V600E, 3: MEK2_C125S, : 4: p61, 5:EGFR, 6:NRAS_Q61K
    ReplicateRows = 2+[(CellLine-1)*4:(CellLine-1)*4+3]; % Choose the appropriate rows for the given cell line
    DoseResponseVec = mean( DoseResponseData(ReplicateRows,2:end),1) ; % Calculate the mean percent suriving across the replications 
    %% Plotting
    Fig1 =  figure(1);
    g60 = scatter(DoseVec,DoseResponseVec,60,ColorMatrix(CellLine, :),'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
    ylim([-0.05 1.05])
    set(gca,'xscale','log')
end





close all
clear all
format long
clc

%% Viral Load data
FittingSwitch =  1; % MEK or BRAF data
if FittingSwitch == 1;
    DoseResponseDataMatrix = struct2array(load('MEKi_DoseResponse_Mutant.mat'));
    DoseResponseData = DoseResponseDataMatrix; % load data
    NReplicates = 4; % how many replicates are there
    CellLine = 1;  % 1 : WT, 2 :V600E, 3: MEK2_C125S, : 4: p61, 5:EGFR, 6:NRAS_Q61K % Which cell line to consider
    ReplicateRows = 2+[(CellLine-1)*4:(CellLine-1)*4+3]; % Choose the appropriate rows for the given cell line
    DoseResponseVec = mean( DoseResponseData(ReplicateRows,2:end),1) ; % Calculate the mean percent suriving across the replications
else
    DoseResponseDataMatrix = struct2array(load('BRAFi_DoseResponse_Mutant.mat'));
    DoseResponseData = DoseResponseDataMatrix; % load data
    NReplicates = 4; % how many replicates are there
    CellLine = 1;  % 1 : WT, 2 :V600E, 3: MEK2_C125S, : 4: p61, 5:EGFR, 6:NRAS_Q61K % Which cell line to consider
    ReplicateRows = 2+[(CellLine-1)*4:(CellLine-1)*4+3]; % Choose the appropriate rows for the given cell line
    DoseResponseVec = mean( DoseResponseData(ReplicateRows,2:end),1) ; % Calculate the mean percent suriving across the replications
end

%% Parameters to be fit
% % PA.rPop = X(1); %Intrinsic growth rate of tumour population [1/day]
% % PA.CarryingCapacity = 1.0e4; % fixed carrying capacity [cells/mL]
% % PA.gammaMax = X(2); % Death rate in maximal drug concentration [1/day]
% % PA.TreatmentHalfEffect = X(3); % [nMol]

% Bounds
lb = [ 0.1,0.05,1e-1];   
ub = [2,5,1e3];
VectorOfInitialGuess = [ 0.35 ,0.65,1e2]'; % initial guess for parameter set

%%  for multi start: 
rng(7)
NumberofStarts = 5;
StartingPoints =  unifrnd(repmat( [0.5.*VectorOfInitialGuess(1:end-1);1.5*VectorOfInitialGuess(end)],1,NumberofStarts),repmat([1.5.*VectorOfInitialGuess(1:end-1);0.5*VectorOfInitialGuess(end)],1,NumberofStarts));

FValVec = zeros(1,NumberofStarts); 
ParameterVec = zeros(length(VectorOfInitialGuess),NumberofStarts); 
tic
    for kk =   1:NumberofStarts; 
        opts = optimset('MaxFunEvals',4000, 'MaxIter',5000,'Display','iter','TolFun', 1e-8,'TolX',1e-8);
        Y0 =  StartingPoints(:,kk) ; %  VectorOfInitialGuess sampled from the multistart vector
        f = @(X)  LogisticExample_Objective(X,DoseResponseVec);  
        [x,fval] = fmincon(f,Y0,[],[],[],[],lb',ub',[],opts) 
        FValVec(kk) = fval;
        ParameterVec(:,kk) = x; 
    end 
    
    t = datetime('now','Format','dMMMyy');
    S = char(t); 
    filename = [S,'_LogisticExample_Fitting_Results_',N]; 

    save(filename)

    toc
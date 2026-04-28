close all
clear all
format long
clc

%% Load data
MutantSwitch = 0;
% Untreated data
UntreatedDataMatrix = struct2array(load('CAT_Monoculture_NoTreatment_Data.mat'));
UntreatedData = UntreatedDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(UntreatedDataMatrix(2+MutantSwitch,1))); % Normalized to fold change

% Docetaxel data
DocetaxelDataMatrix = struct2array(load('CAT_Monoculture_Docetaxel_Data.mat'));
TreatedData = DocetaxelDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(DocetaxelDataMatrix(2+MutantSwitch,1))); % Normalized to fold change
%% Parameters to be fit
% PA.rPopS = X(1); %Intrinsic growth rate of tumour population [1/day]
% PA.rPopR = X(2); %Intrinsic growth rate of sensitive tumour population [1/day]
% % Treatment parameters
% PA.gammaMaxS = X(3); % Death rate in maximal drug concentration [1/day]
% PA.gammaMaxR = X(4);  % Death rate in maximal drug concentration resistant population [1/day]
% PA.rho = X(5);% Initial fraction of each phenotype
% 
% Bounds
lb = [ 0.1,0.1, 0.15, 1e-6, eps];   
ub = [2,2, 2.5 ,2.5, 1-eps];
VectorOfInitialGuess = [0.65, 0.55, 1.05, 0.0001, 0.1]; % initial guess for parameter set

%%  for multi start: 
% rng(7)
NumberofStarts = 1; % 5;
StartingPoints =  [0.65, 0.55, 1.05, 0.0001, 0.1]'; 

FValVec = zeros(1,NumberofStarts); 
ParameterVec = zeros(length(VectorOfInitialGuess),NumberofStarts); 
tic
    for kk =   1:NumberofStarts; 
        opts = optimset('MaxFunEvals',4000, 'MaxIter',5000,'Display','iter','TolFun', 1e-8,'TolX',1e-8);
        Y0 =  StartingPoints(:,kk) ; %  VectorOfInitialGuess sampled from the multistart vector
        f = @(X)  CAT_Example_Objective(X,UntreatedData,TreatedData);  
        [x,fval] = fmincon(f,Y0,[],[],[],[],lb',ub',[],opts) 
        FValVec(kk) = fval;
        ParameterVec(:,kk) = x; 
    end 
    
    t = datetime('now','Format','dMMMyy');
    S = char(t); 
    filename = [S,'_LogisticExample_Fitting_Results_']; 

%     save(filename)

    toc
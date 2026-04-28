close all
clear all
format long
clc
%% Load data
% Untreated data
UntreatedDataMatrix = struct2array(load('CAT_50WT_50M1_Data.mat'));

UntreatedDataCoculture_WT = UntreatedDataMatrix(2,:).*(UntreatedDataMatrix(2,1)./(UntreatedDataMatrix(2,1))); % Normalized to fold change
UntreatedDataCoculture_Mutant = UntreatedDataMatrix(3,:).*(UntreatedDataMatrix(3,1)./(UntreatedDataMatrix(3,1))); % Normalized to fold change
UntreatedTotal = UntreatedDataCoculture_WT + UntreatedDataCoculture_Mutant;

DocetaxelDataCoculture_WT = UntreatedDataMatrix(4,:).*(UntreatedDataMatrix(2,1)./(UntreatedDataMatrix(4,1))); % Normalized to fold change
DocetaxelDataCoculture_Mutant = UntreatedDataMatrix(5,:).*(UntreatedDataMatrix(3,1)./(UntreatedDataMatrix(5,1))); % Normalized to fold change
DocetaxelTotal = DocetaxelDataCoculture_WT + DocetaxelDataCoculture_Mutant;

%% Parameters to be fit
% PA.rPopS = X(1); %Intrinsic growth rate of tumour population [1/day]
% PA.rPopR = X(2); %Intrinsic growth rate of sensitive tumour population [1/day]
% % Treatment parameters
% PA.gammaMaxS = X(3); % Death rate in maximal drug concentration [1/day]
% PA.gammaMaxR = X(4);  % Death rate in maximal drug concentration resistant population [1/day]
% PA.a10 = X(5); % Competition/co-operation  term
% PA.a01 = X(6);  % Competition/co-operation term
% 
% Bounds
lb = [ 0.1,0.1, 0.15, 1e-6, -5,-5];   
ub = [2,2, 2.5 ,2.5, 5, 5];
VectorOfInitialGuess = [0.65, 0.55, 1.05, 0.0001, 1,1]; % initial guess for parameter set

%%  for multi start: 
% rng(7)
NumberofStarts = 1; % 5;
StartingPoints =  VectorOfInitialGuess'; 

FValVec = zeros(1,NumberofStarts); 
ParameterVec = zeros(length(VectorOfInitialGuess),NumberofStarts); 
tic
    for kk =   1:NumberofStarts; 
        opts = optimset('MaxFunEvals',4000, 'MaxIter',5000,'Display','iter','TolFun', 1e-8,'TolX',1e-8);
        Y0 =  StartingPoints(:,kk) ; %  VectorOfInitialGuess sampled from the multistart vector

        f = @(X)  CAT_Example_CoCulture_Objective(X,UntreatedDataCoculture_WT,UntreatedDataCoculture_Mutant,UntreatedTotal,...
                    DocetaxelDataCoculture_WT, DocetaxelDataCoculture_Mutant, DocetaxelTotal);  
        [x,fval] = fmincon(f,Y0,[],[],[],[],lb',ub',[],opts) 
        FValVec(kk) = fval;
        ParameterVec(:,kk) = x; 
    end 
    
    t = datetime('now','Format','dMMMyy');
    S = char(t); 
    filename = [S,'_LogisticExample_Fitting_Results_']; 

%     save(filename)

    toc
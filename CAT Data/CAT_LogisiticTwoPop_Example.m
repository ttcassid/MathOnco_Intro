%% This script simulates the cancer cell dynamics in presence of cytotoxic drug to replicated dose response experiments
close all
clear all
tic
%% Load data
MutantSwitch = 0;
% Untreated data
UntreatedDataMatrix = struct2array(load('CAT_Monoculture_NoTreatment_Data.mat'));
UntreatedData = UntreatedDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(UntreatedDataMatrix(2+MutantSwitch,1))); % Normalized to fold change

% Docetaxel data
DocetaxelDataMatrix = struct2array(load('CAT_Monoculture_Docetaxel_Data.mat'));
TreatedData = DocetaxelDataMatrix(2+MutantSwitch,:).*(UntreatedDataMatrix(2+MutantSwitch,1)./(DocetaxelDataMatrix(2+MutantSwitch,1))); % Normalized to fold change

%% Set model parameters
X = [0.65, 0.55, 1.05, 0.0001, 0.1];
%% Parameters for treatment switching
PA.rPopS = X(1); %Intrinsic growth rate of tumour population [1/day]
PA.rPopR = X(2); %Intrinsic growth rate of sensitive tumour population [1/day]
% Treatment parameters
PA.gammaMaxS = X(3); % Death rate in maximal drug concentration [1/day]
PA.gammaMaxR = X(4);  % Death rate in maximal drug concentration resistant population [1/day]
PA.rho = X(5);% Initial fraction of each phenotype

PA.CarryingCapacity = 1.0e6; % fixed carrying capacity [cells/mL]
PA.TreatmentHalfEffect = 1e2; % [nMol] % Fixed in this setting, normally would be fit [nMol]

%% Treatment Parameters
PA.StartTime = 3;
TimeVec = [0,2,4,6];
NDataPoints = length(TimeVec);
TumourFinalSize = zeros(NDataPoints,1); % Initialize the vector to store final tumour size
Obj = zeros(NDataPoints,1); % Vector to calculate error between model prediction and observed data
%% Time set up
t0 = TimeVec(1); % Start time
tf = TimeVec(end);  %7 day experiment
TotalTime = [t0 tf];
%% Inital Conditions
TumourIC = UntreatedData(1);
SensitiveIC = (1-PA.rho).*TumourIC ; % [cells/mL], Initial condition for the sensitive cells
ResistantIC = (PA.rho).*TumourIC ; % [cells/mL], Initial condition for the sensitive cells
%% Solve the ODE systems for zero dose
TreatCapacityIC = [SensitiveIC,ResistantIC]; % The initial conditions are fixed for all runs
PA.DrugConc = 0;
%% Solve the Untreated ODE system for days 1- 7 
[sol1Untreated] =  FixedDoseDynamicsTreated(TotalTime,TreatCapacityIC,PA); % FixedTherapyDynamicsTreated(TotalTime,TreatCapacityIC,PA); %
TreatedIC =  [deval(sol1Untreated,3,1), deval(sol1Untreated,3,2)];
EvalSol1 = deval(sol1Untreated,TimeVec,1)+ deval(sol1Untreated,TimeVec,2) ; % Calculate population size at each of the measurement times
ObjUntreated = sum( (log10(EvalSol1) - log10(UntreatedData) ).^2 );
%% Solve the treated ODE system
PA.DrugConc = 500 ; 
[sol1Treated] =  FixedDoseDynamicsTreated([3,tf] ,TreatedIC,PA); % FixedTherapyDynamicsTreated(TotalTime,TreatCapacityIC,PA); %
EvalSol2 = deval(sol1Treated,TimeVec([3,4]),1)+ deval(sol1Treated,TimeVec([3,4]),2) ; % Calculate population size at each of the measurement times
ObjTreated = sum( (log10(EvalSol2) - log10(TreatedData([3,4])) ).^2 );

ObjTotal = ObjUntreated + ObjTreated;
toc
 %% Plotting
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

TreatmentColorIndex = 1;

Fig1 =  figure(1); 
hold on
plot(sol1Untreated.x,sol1Untreated.y(1,:)+sol1Untreated.y(2,:),'LineWidth',2,'Color', ColorMatrix(1, :),'LineStyle','-' )
hold on
g60 = scatter(TimeVec,UntreatedData,60,ColorMatrix(1, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
g1 = scatter(TimeVec([3,4]),TreatedData([3,4]),40,ColorMatrix(2, :),'s','LineWidth',1.85);
hold on
plot(sol1Treated.x,sol1Treated.y(1,:)+sol1Treated.y(2,:),'LineWidth',2,'Color', ColorMatrix(2, :),'LineStyle','-' )
hold on
xlim([-0.15 6.15])
% set(gca,'yscale','log')
toc


%% Solvers
function [sol] = FixedDoseDynamicsTreated(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'MaxStep',1e-2);
sol = ode15s(@FixedTherapyDynamics,totaltime,IC,opts);
    function dydt = FixedTherapyDynamics(t,y);
        dydt(1) = PA.rPopS.*y(1).*(1 - (y(1)+y(2))./PA.CarryingCapacity ) -  (PA.gammaMaxS).*y(1).*(PA.DrugConc./(PA.DrugConc+PA.TreatmentHalfEffect) ) ;  %Differential equation for sensitive cells
        dydt(2) = PA.rPopR.*y(2).*(1 - (y(1)+y(2))./PA.CarryingCapacity ) -  (PA.gammaMaxR).*y(2).*(PA.DrugConc./(PA.DrugConc+PA.TreatmentHalfEffect) ) ;  %Differential equation for resistant cells
        dydt = dydt';
    end

end




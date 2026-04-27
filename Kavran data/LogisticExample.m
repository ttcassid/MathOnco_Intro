%% This script simulates the cancer cell dynamics in presence of cytotoxic drug to replicated dose response experiments
close all
clear all
tic
%% Load data
% MEKi data
DoseResponseDataMatrix = struct2array(load('MEKi_DoseResponse_Mutant.mat')); % struct2array(load('BRAFi_DoseResponse_Mutant.mat'));

DoseResponseData = DoseResponseDataMatrix; % load data
NReplicates = 4; % how many replicates are there
CellLine = 1;  % 1 : WT, 2 :V600E, 3: MEK2_C125S, : 4: p61, 5:EGFR, 6:NRAS_Q61K % Which cell line to consider
ReplicateRows = 2+[(CellLine-1)*4:(CellLine-1)*4+3]; % Choose the appropriate rows for the given cell line
DoseResponseVec = mean( DoseResponseData(ReplicateRows,2:end),1) ; % Calculate the mean percent suriving across the replications
%% Set model parameters
X = [0.35, 0.65, 1e2];
%% Parameters for treatment switching
PA.rPop = X(1); %Intrinsic growth rate of tumour population [1/day]
PA.CarryingCapacity = 1.0e4; % fixed carrying capacity [cells/mL]
PA.gammaMax = X(2); % Death rate in maximal drug concentration [1/day]
PA.TreatmentHalfEffect = X(3); % [nMol]
%% Time set up
t0 = 0; % Start time
tf = 3;  %72 hr experiment
TotalTime = [t0 tf];
%% Treatment Parameters
PA.StartTime = t0;
DoseVec = [0.0,0.1,0.3,1,3,10,30,100,300,1000,3000,10000]; % Doses for which there is experimental data
NDataPoints = length(DoseVec);
TumourFinalSize = zeros(NDataPoints,1); % Initialize the vector to store final tumour size
Obj = zeros(NDataPoints,1); % Vector to calculate error between model prediction and observed data
%% Inital Conditions
TumourPopIC = 10; % [cells/mL], Normalized initial condition for dose response curve
%% Solve the ODE systems for zero dose
TreatCapacityIC = [TumourPopIC]; % The initial conditions are fixed for all runs
PA.DrugConc = 0;
%% Solve the Treated ODE system
[sol1Treated] =    FixedDoseDynamicsTreated(TotalTime,TreatCapacityIC,PA); % FixedTherapyDynamicsTreated(TotalTime,TreatCapacityIC,PA); %
TumourFinalSize(1) =  deval(sol1Treated,tf,1);
Obj(1) = ( log10(1) - log10(DoseResponseVec(1)) ).^2;

%% Solve the ODE systems for non-zero dose
for ii = 2:  NDataPoints
    PA.DrugConc = DoseVec(ii);
    %% Solve the Treated ODE system
    [sol1Treated] =    FixedDoseDynamicsTreated(TotalTime,TreatCapacityIC,PA); % FixedTherapyDynamicsTreated(TotalTime,TreatCapacityIC,PA); %
    TumourFinalSize(ii) =  deval(sol1Treated,tf,1)./TumourFinalSize(1); % Relative growth compared to no dosing
    
    Obj(ii) = ( log10(TumourFinalSize(ii)) - log10(DoseResponseVec(ii)) ).^2; % Error on log10 scale between model prediction and outcome
end
ObjTotal = sum(Obj)
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
g60 = scatter(DoseVec(2:end),DoseResponseVec(2:end),60,ColorMatrix(CellLine, :),'o','filled'); %,'20', [33,102,172]/255,'*');
hold on
g1 = scatter(DoseVec(2:end),TumourFinalSize(2:end),40,ColorMatrix(10, :),'s','LineWidth',1.85);
hold on
ylim([-0.05 1.05])
set(gca,'xscale','log')

toc


%% Solvers
function [sol] = FixedDoseDynamicsTreated(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'MaxStep',1e-2);
sol = ode15s(@FixedTherapyDynamics,totaltime,IC,opts);
    function dydt = FixedTherapyDynamics(t,y);
        dydt(1) = PA.rPop.*y(1).*(1 - (y(1))./PA.CarryingCapacity ) -  (PA.gammaMax).*(PA.DrugConc./(PA.DrugConc+PA.TreatmentHalfEffect) ).*y(1) ; ; %Differential equation for N(t)
        dydt = dydt';
    end

end

function Dose = DoseTherapy(PA,t); %IV administration of therapy
DoseVec = zeros(1,PA.AdminNumber); %Create a vector of doses for each administration
for nn = 1:PA.AdminNumber
    if PA.TreatmentStartTime(nn)  < t && t < PA.TreatmentEndTime(nn)
        DoseVec(nn) = (PA.Admin(nn))./(PA.Vol.*PA.TimeAdmin);
    else
        DoseVec(nn) = 0;
    end
end
Dose = ones(1,PA.AdminNumber)*DoseVec';
end



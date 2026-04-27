function [ObjTotal] = LogisticExample_Objective(X,DoseResponseVec);
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
DoseVec = [0.01,0.1,0.3,1,3,10,30,100,300,1000,3000,10000]; % Doses for which there is experimental data
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
Obj(1) = ( log10(TumourFinalSize(1)) - log10(DoseResponseVec(1)) ).^2;

%% Solve the ODE systems for non-zero dose
for ii = 2:  NDataPoints
    PA.DrugConc = DoseVec(ii);
    %% Solve the Treated ODE system
    [sol1Treated] =    FixedDoseDynamicsTreated(TotalTime,TreatCapacityIC,PA); % FixedTherapyDynamicsTreated(TotalTime,TreatCapacityIC,PA); %
    TumourFinalSize(ii) =  deval(sol1Treated,tf,1)./TumourFinalSize(1); % Relative growth compared to no dosing
    
    Obj(ii) = ( log10(TumourFinalSize(ii)) - log10(DoseResponseVec(ii)) ).^2; % Error on log10 scale between model prediction and outcome
end 
ObjTotal = sum(Obj);
end

function [sol] = FixedDoseDynamicsTreated(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'MaxStep',1e-2);
sol = ode15s(@FixedTherapyDynamics,totaltime,IC,opts);
    function dydt = FixedTherapyDynamics(t,y);
        dydt(1) = PA.rPop.*(1 - (y(1))./PA.CarryingCapacity ) -  (PA.gammaMax).*(PA.DrugConc./(PA.DrugConc+PA.TreatmentHalfEffect) ).*y(1) ; ; %Differential equation for N(t)
        dydt = dydt';
    end

end
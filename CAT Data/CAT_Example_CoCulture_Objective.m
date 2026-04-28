function [ObjTotal] = CAT_Example_CoCulture_Objective(X,UntreatedDataCoculture_WT,UntreatedDataCoculture_Mutant,UntreatedTotal,...
                    DocetaxelDataCoculture_WT, DocetaxelDataCoculture_Mutant, DocetaxelTotal);
%% Parameters for treatment switching
PA.rPopS = X(1); %Intrinsic growth rate of tumour population [1/day]
PA.rPopR = X(2); %Intrinsic growth rate of sensitive tumour population [1/day]
% Treatment parameters
PA.gammaMaxS = X(3); % Death rate in maximal drug concentration [1/day]
PA.gammaMaxR = X(4);  % Death rate in maximal drug concentration resistant population [1/day]
PA.a10 = X(5); % Competition/co-operation  term
PA.a01 = X(6);  % Competition/co-operation term


PA.CarryingCapacity = 1.0e6; % fixed carrying capacity [cells/mL]
PA.TreatmentHalfEffect = 1e2; % [nMol] % Fixed in this setting, normally would be fit [nMol]
PA.rho = 0.5;% Initial fraction of each Genotype
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
TumourIC = UntreatedTotal(1);
WTIC = (1-PA.rho).*TumourIC ; % [cells/mL], Initial condition for the sensitive cells
M1IC = (PA.rho).*TumourIC ; % [cells/mL], Initial condition for the sensitive cells
%% Solve the ODE systems for zero dose
TreatCapacityIC = [WTIC,M1IC]; % The initial conditions are fixed for all runs
PA.DrugConc = 0;
%% Solve the Untreated ODE system for days 1- 7 
[sol1Untreated] =  FixedDoseDynamicsTreated(TotalTime,TreatCapacityIC,PA);  
TreatedIC =  [deval(sol1Untreated,3,1), deval(sol1Untreated,3,2)];
EvalSol1 = deval(sol1Untreated,TimeVec,1)+ deval(sol1Untreated,TimeVec,2) ; % Calculate population size at each of the measurement times
ObjUntreated = sum( (log10(EvalSol1) - log10(UntreatedTotal) ).^2 ) ...
    + sum(log10(deval(sol1Untreated,TimeVec,1)) - log10(UntreatedDataCoculture_WT) ).^2 ...
    + sum(log10(deval(sol1Untreated,TimeVec,2)) - log10(UntreatedDataCoculture_Mutant) ).^2  ;

%% Solve the treated ODE system
PA.DrugConc = 500 ; 
[sol1Treated] =  FixedDoseDynamicsTreated([3,tf] ,TreatedIC,PA); % FixedTherapyDynamicsTreated(TotalTime,TreatCapacityIC,PA); %
EvalSol2 = deval(sol1Treated,TimeVec([3,4]),1)+ deval(sol1Treated,TimeVec([3,4]),2) ; % Calculate population size at each of the measurement times
ObjTreated = sum( (log10(EvalSol2) - log10(DocetaxelTotal([3,4])) ).^2 ) ...
    + sum(log10(deval(sol1Treated,TimeVec([3,4]),1)) - log10(DocetaxelDataCoculture_WT([3,4])) ).^2 ... 
    + sum(log10(deval(sol1Treated,TimeVec([3,4]),2)) - log10(DocetaxelDataCoculture_Mutant([3,4])) ).^2  ;

ObjTotal = ObjUntreated + ObjTreated;

end

%% Solvers
function [sol] = FixedDoseDynamicsTreated(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-7,'AbsTol',1e-7,'MaxStep',1e-2);
sol = ode15s(@FixedTherapyDynamics,totaltime,IC,opts);
    function dydt = FixedTherapyDynamics(t,y);
        dydt(1) = PA.rPopS.*y(1).*(1 - (y(1)+PA.a10.*y(2))./PA.CarryingCapacity ) -  (PA.gammaMaxS).*y(1).*(PA.DrugConc./(PA.DrugConc+PA.TreatmentHalfEffect) ) ;  %Differential equation for sensitive cells
        dydt(2) = PA.rPopR.*y(2).*(1 - (y(1)+PA.a01.*y(2))./PA.CarryingCapacity ) -  (PA.gammaMaxR).*y(2).*(PA.DrugConc./(PA.DrugConc+PA.TreatmentHalfEffect) ) ;  %Differential equation for resistant cells
        dydt = dydt';
    end

end
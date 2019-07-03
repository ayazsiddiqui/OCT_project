clear
clc
close all
format compact

plot_animation = 0;
make_video = 0;

%%
sim_time = 500;
tVec = 0:0.05:sim_time;

% setpoints

altitudeSP = 100*ones(size(tVec));

altiMax = 175;
altiMin = 125;

pitchSP = 8*(pi/180)*ones(size(tVec));

rollPeriod = 100;

rollAmp = 20;
rollSP = (pi/180)*rollAmp*sign(sin((2*pi/rollPeriod)*tVec));

altitudeSP = timeseries(altitudeSP,tVec);
pitchSP = timeseries(pitchSP,tVec);
rollSP = timeseries(rollSP,tVec);


%% environment conditions
env_t = sysParam.env;

env_t.gravAccel.value = 9.81;
env_t.flowDensity.value = 1000;
env_t.inertialFlowVel.value = [1;0;0];

% created new branch

% create class instance
tp = sysParam.plant_v2(1,2);
tp.lengthScaleFactor = 1;
tp.densityScaleFactor = 1;
tp.buoyancyFactor = 1.1;

%% set vehicle values

tp.vehicle.MI.value = 1e-6*[6.303080401918E+09 0 0;...
    0 2080666338.077 0;...
    0 0 8.320369733598E+09];
tp.vehicle.volume.value = 945352023e-9;
tp.vehicle.mass.value = (1/tp.buoyancyFactor)*tp.vehicle.volume.value*env_t.flowDensity.value;

tp.vehicle.Rcb_cm.value = [0;0;0];
tp.vehicle.Rcm_wingLE.value = [1.5;0;0];

% initial operating conditions
tp.vehicle.ini_Rcm_o.value = [0; 0; altitudeSP.Data(1)];
tp.vehicle.ini_O_Vcm_o.value = env_t.inertialFlowVel.value;
tp.vehicle.ini_euler.value = [0; 4; 0]*pi/180;
tp.vehicle.ini_OwB.value = [0; 0; 0];

%% partitioned lifting body parameters
tp.aeroDataFileName = 'partDsgn1_lookupTables.mat';

tp = tp.calcAddedMass(env_t);

%% turbines
tp.turbines(1).Rturb_cm.value = [2.5000; -20.2500; 0];
tp.turbines(1).diameter.value = 0;
tp.turbines(1).powerCoeff.value = 0.5;
tp.turbines(1).dragCoeff.value = 0.8;

tp.turbines(2).Rturb_cm.value = [2.5000; 20.2500; 0];
tp.turbines(2).diameter.value = 0;
tp.turbines(2).powerCoeff.value = 0.5;
tp.turbines(2).dragCoeff.value = 0.8;

%% tethers
tp.tethers(1).numNodes = 4;
tp.tethers(1).diameter = 0.01;
tp.tethers(1).youngsModulus = 3.8e9;
tp.tethers(1).dampingRatio = 0.05;
tp.tethers(1).dragCoeff = 0.5;
tp.tethers(1).density = 1300;

% redesign tethers
maxAppFlowMultiplier = 4;
maxPercentageElongation = 0.02;

tp = tp.designTetherDiameter(env_t,maxAppFlowMultiplier,maxPercentageElongation);

%% gnd station
% rotation switch
tp.gndStation.rotationSwitch.value = 0;

% parameters
tp.gndStation.Izz.value = 100;
tp.gndStation.dampCoeff.value = 10;

% initial conditions
tp.gndStation.ini_platform_ang.value = 0*pi/180;
tp.gndStation.ini_platform_vel.value = 0;

%% winches
tp.winches(1).maxSpeed = 1;
tp.winches(1).timeConstant = 0.1;


tp = tp.setTetherInitLength(env_t);

%% define controller parameter
% tether command gains
ctrllr.tethers.transformMat = [1 .5 -.5; 1 -.5 0; 1 .5 .5];


ctrllr.tethers.altiTetherKp = 0*0.5;    % m/s per m
ctrllr.tethers.altiTetherKi = 0;
ctrllr.tethers.altiTetherKd = 0*0.1;
ctrllr.tethers.altiTetherTau = 1;
ctrllr.tethers.altiErrorSat = 1;


ctrllr.tethers.pitchTetherKp = 0.1;   % m/s per rad
ctrllr.tethers.pitchTetherKi = 0;
ctrllr.tethers.pitchTetherKd = 0.2;
ctrllr.tethers.pitchTetherTau = 0.8;

ctrllr.tethers.rollTetherKp = 4;    % m/s per rad
ctrllr.tethers.rollTetherKi = 0;
ctrllr.tethers.rollTetherKd = 12;      % m/s per rad/s
ctrllr.tethers.rollTetherTau = 0.5;

% control surface gains
ctrllr.controlSurfaces.aileronKp = 1*0.75;   % deg per deg
ctrllr.controlSurfaces.aileronKi = 0;
ctrllr.controlSurfaces.aileronKd = 1.5;
ctrllr.controlSurfaces.aileronTau = 0.5;
ctrllr.controlSurfaces.aileronMaxDef = 30;

ctrllr.controlSurfaces.elevatorKp = 1*1;  % deg per deg
ctrllr.controlSurfaces.elevatorKi = 0.00;
ctrllr.controlSurfaces.elevatorKd = 0.2;
ctrllr.controlSurfaces.elevatorTau = 0.5;
ctrllr.controlSurfaces.elevatorMaxDef = 30;

%% simulate
run_no = 1;

% max step
maxStep = 0.1;
simWithMonitor('mainModel',2)
save('unscaled_res')
scaledModel_postProcess

if plot_animation == 1
    fullKitePlotter
end

% scaleFactors(1) = 1/100;
% scaleFactors(2) = 1;
% 
% [s_tp,s_env_t,s_ctrllr,s_sim_time,s_altitudeSP,s_pitchSP,s_rollSP] = ...
%     scaleEverything(scaleFactors,tp,env_t,ctrllr,sim_time,altitudeSP,pitchSP,rollSP);
% 
% clear tp env_t ctrllr sim_time altitudeSP pitchSP rollSP
% 
% tp = s_tp;
% env_t = s_env_t;
% ctrllr = s_ctrllr;
% sim_time = s_sim_time;
% altitudeSP = s_altitudeSP;
% pitchSP = s_pitchSP;
% rollSP = s_rollSP;
% 
% run_no = run_no + 1;
% 
% simWithMonitor('mainModel',2)
% save('scaled_res')
% scaledModel_postProcess

clear
clc
format compact

%%
%% environment conditions
env_t = sysParam.env;

env_t.gravAccel.value = 9.81;
env_t.flowDensity.value = 1000;
env_t.iniertialFlowVel.value = [1;0;0];

% created new branch

% create class instance
tp = sysParam.plant_v2(3,2);
tp.ScaleFactor = 1;

%% set vehicle values
tp.vehicle.mass.value = 945.352/1.05;
tp.vehicle.MI.value = 1e-6*[6.303080401918E+09 0 0;...
    0 2080666338.077 0;...
    0 0 8.320369733598E+09];
tp.vehicle.volume.value = 945352023e-9;
tp.vehicle.Rcb_cm.value = [0;0;0];
tp.vehicle.Rcm_wingLE.value = [1.5;0;0];

% setpoints
altitudeSP = 200;
pitchSP = 2*pi/180;
rollSP = 0*pi/180;

% initial operating conditions
tp.vehicle.ini_Rcm_o.value = [0; 0; altitudeSP];
tp.vehicle.ini_O_Vcm_o.value = [0; 0; 0];
tp.vehicle.ini_euler.value = [0; pitchSP; 0];
tp.vehicle.ini_OwB.value = [0; 0; 0];

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
tp.tethers(1).diameter = 0.055;
tp.tethers(1).youngsModulus = 3.8e9;
tp.tethers(1).dampingRatio = 0.05;
tp.tethers(1).dragCoeff = 0.5;
tp.tethers(1).density = 1300;

tp.tethers(2).numNodes = 4;
tp.tethers(2).diameter = 0.076;
tp.tethers(2).youngsModulus = 3.8e9;
tp.tethers(2).dampingRatio = 0.05;
tp.tethers(2).dragCoeff = 0.5;
tp.tethers(2).density = 1300;

tp.tethers(3).numNodes = 4;
tp.tethers(3).diameter = 0.055;
tp.tethers(3).youngsModulus = 3.8e9;
tp.tethers(3).dampingRatio = 0.05;
tp.tethers(3).dragCoeff = 0.5;
tp.tethers(3).density = 1300;

%% partitioned lifting body parameters
tp.aeroDataFileName = 'partDsgn1_lookupTables.mat';

tp = tp.calcAddedMass(env_t);

% redesign tethers
maxAppFlowMultiplier = 2;
maxPercentageElongation = 0.05;

tp = tp.designTetherDiameter(env_t,maxAppFlowMultiplier,maxPercentageElongation);

%% winches
tp.winches(1).maxSpeed = 0.4;
tp.winches(1).timeConstant = 1;
tp.winches(1).initTetherLength = 199.0076;

tp.winches(2).maxSpeed = 0.4;
tp.winches(2).timeConstant = 1;
tp.winches(2).initTetherLength = 196.0955;

tp.winches(3).maxSpeed = 0.4;
tp.winches(3).timeConstant = 1;
tp.winches(3).initTetherLength = 199.0076;

%% gnd station
tp.gndStation.Izz.value = 100;
tp.gndStation.dampCoeff.value = 10;


% initial conditions
tp.gndStation.ini_platform_ang.value = 0*pi/180;
tp.gndStation.ini_platform_vel.value = 0;

tp = tp.setTetherInitLength(env_t);


%% define controller parameter
% tether command gains
ctrllr.tethers.transformMat = [1 .5 -.5; 1 -.5 0; 1 .5 .5];

ctrllr.tethers.altiTetherKp = 0.05;    % m/s per m
ctrllr.tethers.altiTetherKi = 0;
ctrllr.tethers.altiTetherKd = 5*ctrllr.tethers.altiTetherKp;
ctrllr.tethers.altiTetherTau = 31.8310;

ctrllr.tethers.pitchTetherKp = 1.5*0.5;   % m/s per rad
ctrllr.tethers.pitchTetherKi = 0;
ctrllr.tethers.pitchTetherKd = 2.5*ctrllr.tethers.pitchTetherKp;
ctrllr.tethers.pitchTetherTau = 0.7958;

ctrllr.tethers.rollTetherKp = 1*ctrllr.tethers.pitchTetherKp/(0.5*tp.aeroDesignData.wing_AR);    % m/s per rad
ctrllr.tethers.rollTetherKi = 0;
ctrllr.tethers.rollTetherKd = 2*ctrllr.tethers.rollTetherKp;
ctrllr.tethers.rollTetherTau = 0.7958;

% control surface gains
ctrllr.controlSurfaces.aileronKp = 2;   % deg per deg
ctrllr.controlSurfaces.aileronKi = 0;
ctrllr.controlSurfaces.aileronKd = 2*ctrllr.controlSurfaces.aileronKp;
ctrllr.controlSurfaces.aileronTau = 0.2;
ctrllr.controlSurfaces.aileronMaxDef = 30;

ctrllr.controlSurfaces.elevatorKp = 1;  % deg per deg
ctrllr.controlSurfaces.elevatorKi = 0;
ctrllr.controlSurfaces.elevatorKd = 3*ctrllr.controlSurfaces.elevatorKp;
ctrllr.controlSurfaces.elevatorTau = 0.05;

ctrllr.controlSurfaces.elevatorMaxDef = 30;










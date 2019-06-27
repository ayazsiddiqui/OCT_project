clear
clc
format compact

%% environment conditions
flowVel = [0.2;0;0];
rho_fluid = 1000;
grav = 9.81;

% create class instance
tp = plant(3,2);
tp.ScaleFactor = 1;

%% set vehicle values
tp.vehicle.mass.value = 8.9360e+04;
tp.vehicle.MI.value = [14330000 0 0; 0 143200 0; 0 0 15300000];
tp.vehicle.volume.value = 111.7000;
tp.vehicle.Rcb_cm.value = [0;0;0];
tp.vehicle.Rcm_wingLE.value = [0;0;0];

% attachment points
tp.vehicle.tetherAttchPts(1).value = [-2.5000; -20.0000; -0.3750];
tp.vehicle.tetherAttchPts(2).value = [21.2500; 0; -0.3750];
tp.vehicle.tetherAttchPts(3).value = [-2.5000; 20.0000; -0.3750];

% setpoints
altitudeSP = 200;
pitchSP = 2*pi/180;
rollSP = 0*pi/180;

% initial operating conditions
ini_Rcm_o = [0; 0; 200];
ini_O_Vcm_o = [0; 0; 0];
ini_euler = [0; 0.1222; 0];
ini_OwB = [0; 0; 0];


%% turbines
tp.turbines(1).Rturb_cm.value = [2.5000; -20.2500; 0];
tp.turbines(1).diameter.value = 8.7;
tp.turbines(1).powerCoeff.value = 0.5;
tp.turbines(1).dragCoeff.value = 0.8;

tp.turbines(2).Rturb_cm.value = [2.5000; 20.2500; 0];
tp.turbines(2).diameter.value = 8.7;
tp.turbines(2).powerCoeff.value = 0.5;
tp.turbines(2).dragCoeff.value = 0.8;

%% tethers
tp.tethers(1).numNodes.value = 4;
tp.tethers(1).diameter.value = 0.055;
tp.tethers(1).youngsModulus.value = 3.8e9;
tp.tethers(1).dampingRation.value = 0.05;
tp.tethers(1).dragCoeff.value = 0.5;
tp.tethers(1).density.value = 1300;
tp.tethers(1).vehicleMass.value = tp.vehicle.mass.value;

tp.tethers(2).numNodes.value = 4;
tp.tethers(2).diameter.value = 0.076;
tp.tethers(2).youngsModulus.value = 3.8e9;
tp.tethers(2).dampingRation.value = 0.05;
tp.tethers(2).dragCoeff.value = 0.5;
tp.tethers(2).density.value = 1300;
tp.tethers(2).vehicleMass.value = tp.vehicle.mass.value;

tp.tethers(3).numNodes.value = 4;
tp.tethers(3).diameter.value = 0.055;
tp.tethers(3).youngsModulus.value = 3.8e9;
tp.tethers(3).dampingRation.value = 0.05;
tp.tethers(3).dragCoeff.value = 0.5;
tp.tethers(3).density.value = 1300;
tp.tethers(3).vehicleMass.value = tp.vehicle.mass.value;

%% winches
tp.winches(1).maxSpeed.value = 0.4;
tp.winches(1).timeConstant.value = 1;
tp.winches(1).initTetherLength.value = 199.0076;

tp.winches(2).maxSpeed.value = 0.4;
tp.winches(2).timeConstant.value = 1;
tp.winches(2).initTetherLength.value = 196.0955;

tp.winches(3).maxSpeed.value = 0.4;
tp.winches(3).timeConstant.value = 1;
tp.winches(3).initTetherLength.value = 199.0076;

%% gnd station
tp.gndStation.Izz.value = 100;
tp.gndStation.dampCoeff.value = 10;
tp.gndStation.tetherAttchPts(1).value = [-2.5000; -20.0000; 0.3750];
tp.gndStation.tetherAttchPts(2).value = [21.2500; 0; 0.3750];
tp.gndStation.tetherAttchPts(3).value = [21.2500; 0; 0.3750];

% initial conditions
ini_platform_ang = 0;
ini_platform_vel = 0;

%% partitioned lifting body parameters
load('partDsgn1_lookupTables.mat')

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

ctrllr.tethers.rollTetherKp = 1*ctrllr.tethers.pitchTetherKp/(0.5*dsgnData.wing_AR);    % m/s per rad
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















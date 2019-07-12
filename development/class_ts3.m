clear
clc
format compact
% close all

% % merged

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

%% environment
env = ENV.environment;

env.setLengthScale(lengthScale,'');
env.setDensityScale(densityScale,'');

env.setInertialFlowVel([1;0;0],'m/s');

env.scaleEnvironment;

%% lifiting body
vhcl = PLT.vehicle;

vhcl.setLengthScale(lengthScale,'');
vhcl.setDensityScale(densityScale,'');
vhcl.setNumTethers(numTethers,'');
vhcl.setNumTurbines(numTurbines,'');
vhcl.setBuoyFactor(1.25,'');

% % % volume and inertias
vhcl.setVolume(945352023.474*1e-9,'m^3');
vhcl.setIxx(6.303080401918E+09*1e-6,'kg*m^2');
vhcl.setIyy(2080666338.077*1e-6,'kg*m^2');
vhcl.setIzz(8.320369733598E+09*1e-6,'kg*m^2');
vhcl.setIxy(0,'kg*m^2');
vhcl.setIxz(81875397.942*1e-6,'kg*m^2');
vhcl.setIyz(0,'kg*m^2');
vhcl.setRcb_cm([0;0;0],'m');

% % % data file name
vhcl.setFluidCoeffsFileName('somefile','');

% % % wing
vhcl.setRwingLE_cm([-1.0;0;0],'m');
vhcl.setWingChord(1,'m');
vhcl.setWingAR(10,'');
vhcl.setWingTR(0.8,'');
vhcl.setWingSweep(10,'deg');
vhcl.setWingDihedral(2,'deg');
vhcl.setWingIncidence(0,'deg');
vhcl.setWingNACA('2412','');
vhcl.setWingClMax(1.7,'');
vhcl.setWingClMin(-1.7,'');

% % % H-stab
vhcl.setRhsLE_wingLE([6;0;0],'m');
vhcl.setHsChord(0.5,'m');
vhcl.setHsAR(8,'');
vhcl.setHsTR(0.8,'');
vhcl.setHsSweep(10,'deg');
vhcl.setHsDihedral(0,'deg');
vhcl.setHsIncidence(0,'deg');
vhcl.setHsNACA('0015','');
vhcl.setHsClMaxl(1.7,'');
vhcl.setHsClMin(-1.7,'');

% % % V-stab
vhcl.setRvs_wingLE([6;0;0],'m');
vhcl.setVsChord(0.5,'m');
vhcl.setVsSpan(2.5,'m');
vhcl.setVsTR(0.8,'');
vhcl.setVsSweep(10,'deg');
vhcl.setVsNACA('0015','');
vhcl.setVsClMax(1.7,'');
vhcl.setVsClMin(-1.7,'');

% % % initial conditions
vhcl.setInitialCmPos([0;0;100],'m');
vhcl.setInitialCmVel([0;0;0],'m/s');
vhcl.setInitialEuler([0;0;0],'rad');
vhcl.setInitialAngVel([0;0;0],'rad/s');

% % % scale the vehicle
vhcl.scaleVehicle

% % % load/generate fluid dynamic data
vhcl.fluidDynamicCoefffs

% % % plot
% vhcl.plot
% vhcl.plotCoeffPolars

%% turbines
turb = PLT.turbine;

turb.setLengthScale(lengthScale,'');
turb.setDensityScale(densityScale,'');
turb.setNumTurbines(numTurbines,'');
turb.setTurbDiameter([10 10],'m')
turb.setTurbDragCoeff([0.8 0.8],'');
turb.setTurbPowerCoeff([0.5 0.5],'');

% % % scale turbine
turb.scaleTurbine

%% ground station
gnd = PLT.gndStn;

gnd.setLengthScale(lengthScale,'');
gnd.setDensityScale(densityScale,'');
gnd.setNumTethers(numTethers,'');

gnd.setIzz(100,'kg*m^2');
gnd.setDampingCoeff(10,'N*m*s');
gnd.setFreeSpinSwitch(0,'');

gnd.setThrAttchPts(vhcl);

% % % initial conditions
gnd.setInitialEuler(0,'rad');
gnd.setInitialAngVel(0,'rad/s');

gnd.scaleGndStn;

%% tethers
thr = PLT.tether;

thr.setLengthScale(lengthScale,'');
thr.setDensityScale(densityScale,'');
thr.setNumTethers(numTethers,'');

thr.setNumNodes(4,'');
thr.setThrDiameter([0.01 0.02 0.01],'m');
thr.setThrDensity(1300*ones(1,numTethers),'kg/m^3');
thr.setThrYoungs(4e9*ones(1,numTethers),'N/m^2');
thr.setThrDampingRatio(0.05*ones(1,numTethers),'');
thr.setThrDragCoeff(0.5*ones(1,numTethers),'');

thr.scaleTether;

%% winches
wnch = PLT.winch;

wnch.setLengthScale(lengthScale,'');
wnch.setDensityScale(densityScale,'');
wnch.setNumTethers(numTethers,'');

wnch.setWnchMaxTugSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchMaxReleaseSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchTimeConstant(1*ones(1,numTethers),'s');

wnch.scaleWinch;




















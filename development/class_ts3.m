clear
clc
format compact
% close all

% % merged

%% common parameters
lengthScale = 0.5;
densityScale = 1;
numTethers = 3;
numTurbines = 2;


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
turb.setTurbDiameter([0 0],'m')
turb.setTurbDragCoeff([0.8 0.8],'');
turb.setTurbPowerCoeff([0.5 0.5],'');

% % % scale turbine
turb.scaleTurbine

%% ground station












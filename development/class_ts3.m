clear
clc
format compact
% close all

%% simtime
plot_animation = 0;
make_video = 0;

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

sim_time = 200*sqrt(lengthScale);

%% set variants
vhcl_variant = 'partitionedLiftingBodyVariant';
thr_variant = 'KelvinVoigtTetherVariant';
wnch_variant = 'PureIntegratorWinchVariant';
gnd_variant = 'FixedGrounStationVariant';

dts = 0.05*sqrt(lengthScale);
tVec = 0:dts:sim_time;

%% setpoints

% altitude
altitudeSP = 100*ones(size(tVec)).*lengthScale;

% pitch
pitchSP = 2*(pi/180)*ones(size(tVec));

% roll
rollAmp = 0;
rollPeriod = 120;
rollSP = (pi/180)*rollAmp*sign(sin((2*pi/rollPeriod)*tVec));

% yaw
yawSP = 0*(pi/180)*ones(size(tVec));

altitudeSP = timeseries(altitudeSP,tVec);
pitchSP = timeseries(pitchSP,tVec);
rollSP = timeseries(rollSP,tVec);
yawSP = timeseries(yawSP,tVec);


%% environment
env = ENV.environment;

env.setLengthScale(lengthScale,'');
env.setDensityScale(densityScale,'');

env.setInertialFlowVel([1 0 0],'m/s');

env.scaleEnvironment;

%% lifiting body
vhcl = PLT.vehicle;

vhcl.setLengthScale(lengthScale,'');
vhcl.setDensityScale(densityScale,'');
vhcl.setNumTethers(numTethers,'');
vhcl.setNumTurbines(numTurbines,'');
vhcl.setBuoyFactor(1.2,'');

% % % volume and inertias
vhcl.setVolume(945352023.474*1e-9,'m^3');
vhcl.setIxx(6.303080401918E+09*1e-6,'kg*m^2');
vhcl.setIyy(2080666338.077*1e-6,'kg*m^2');
vhcl.setIzz(8.320369733598E+09*1e-6,'kg*m^2');
vhcl.setIxy(0,'kg*m^2');
vhcl.setIxz(81875397.942*1e-6,'kg*m^2');
vhcl.setIyz(0,'kg*m^2');
vhcl.setRcb_cm([0;0;0],'m');

% % % wing
vhcl.setRwingLE_cm([-1.0;0;0],'m');
vhcl.setWingChord(1,'m');
vhcl.setWingAR(10,'');
vhcl.setWingTR(0.8,'');
vhcl.setWingSweep(10,'deg');
vhcl.setWingDihedral(2,'deg');
vhcl.setWingIncidence(0,'deg');
vhcl.setWingNACA('0015','');
vhcl.setWingClMax(1.25,'');
vhcl.setWingClMin(-1.25,'');

% % % H-stab
vhcl.setRhsLE_wingLE([6;0;0],'m');
vhcl.setHsChord(0.5,'m');
vhcl.setHsAR(8,'');
vhcl.setHsTR(0.8,'');
vhcl.setHsSweep(10,'deg');
vhcl.setHsDihedral(0,'deg');
vhcl.setHsIncidence(0,'deg');
vhcl.setHsNACA('0015','');
vhcl.setHsClMaxl(1.25,'');
vhcl.setHsClMin(-1.25,'');

% % % V-stab
vhcl.setRvs_wingLE([6;0;0],'m');
vhcl.setVsChord(0.5,'m');
vhcl.setVsSpan(2.5,'m');
vhcl.setVsTR(0.8,'');
vhcl.setVsSweep(10,'deg');
vhcl.setVsNACA('0015','');
vhcl.setVsClMax(1.25,'');
vhcl.setVsClMin(-1.25,'');

% % % initial conditions
vhcl.setInitialCmPos([0;0;100],'m');
vhcl.setInitialCmVel([0;0;0],'m/s');
vhcl.setInitialEuler([0;1;0]*pi/180,'rad');
vhcl.setInitialAngVel([0;0;0],'rad/s');

% % % scale the vehicle
vhcl.scaleVehicle

% % % data file name
vhcl.setFluidCoeffsFileName('somefile1','');

% % % load/generate fluid dynamic data
vhcl.calcFluidDynamicCoefffs

% % % plot
% vhcl.plot
% vhcl.plotCoeffPolars

%% turbines
turb = PLT.turbine;

turb.setLengthScale(lengthScale,'');
turb.setDensityScale(densityScale,'');
turb.setNumTurbines(numTurbines,'');
turb.setTurbDiameter(0*ones(1,numTurbines),'m')
turb.setTurbDragCoeff(0.8*ones(1,numTurbines),'');
turb.setTurbPowerCoeff(0.5*ones(1,numTurbines),'');

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
thr.setThrDiameter(0.01*[1 sqrt(2) 1],'m');
thr.setThrDensity(1300*ones(1,numTethers),'kg/m^3');
thr.setThrYoungs(3.8e9*ones(1,numTethers),'N/m^2');
thr.setThrDampingRatio(0.02*ones(1,numTethers),'');
thr.setThrDragCoeff(0.5*ones(1,numTethers),'');

thr.scaleTether;

% design tether
maxAppFlowMultiplier = 2;
maxPercentageElongation = 0.05;
thr.designTetherDiameter(vhcl,env,maxAppFlowMultiplier,maxPercentageElongation);

%% winches
wnch = PLT.winch;

wnch.setLengthScale(lengthScale,'');
wnch.setDensityScale(densityScale,'');
wnch.setNumTethers(numTethers,'');

wnch.setWnchMaxTugSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchMaxReleaseSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchTimeConstant(1*ones(1,numTethers),'s');

wnch.calcInitTetherLength(vhcl,gnd,env);

wnch.scaleWinch;

%% controller
ctrl = CTR.threeThrCtlr;

ctrl.setLengthScale(lengthScale,'');
ctrl.setDensityScale(densityScale,'');
ctrl.setNumTethers(numTethers,'');

% altitude tether control gains
ctrl.setAltiTetherKp(0.0,'(m/s)/m')
ctrl.setAltiTetherKi(0,'(m/s)/(m*s)')
ctrl.setAltiTetherKd(0,'(m/s)/(m/s)')
ctrl.setAltiTetherTau(1,'s')
ctrl.setAltiErrorSat(5,'m')

% pitch tether control gains
ctrl.setPitchTetherKp(0.1,'(m/s)/rad')
ctrl.setPitchTetherKi(0,'(m/s)/(rad*s)')
ctrl.setPitchTetherKd(0.0,'(m/s)/(rad/s)')
ctrl.setPitchTetherTau(1,'s')

% roll tether control gains
ctrl.setRollTetherKp(0.0,'(m/s)/rad')
ctrl.setRollTetherKi(0,'(m/s)/(rad*s)')
ctrl.setRollTetherKd(0,'(m/s)/(rad/s)')
ctrl.setRollTetherTau(1,'s')

% aileron gains
ctrl.setAileronKp(0,'deg/deg');
ctrl.setAileronKi(0,'deg/(deg*s)');
ctrl.setAileronKd(0,'deg/(deg/s)');
ctrl.setAileronTau(0.5,'s');
ctrl.setAileronMaxDef(30,'deg');

% elevator gains
ctrl.setElevatorKp(0,'deg/deg');
ctrl.setElevatorKi(0,'deg/(deg*s)');
ctrl.setElevatorKd(0,'deg/(deg/s)');
ctrl.setElevatorTau(0.5,'s');
ctrl.setElevatorMaxDef(30,'deg');

% rudder gains
ctrl.setRudderKp(0,'deg/deg');
ctrl.setRudderKi(0,'deg/(deg*s)');
ctrl.setRudderKd(0,'deg/(deg/s)');
ctrl.setRudderTau(0.5,'s');
ctrl.setRudderMaxDef(30,'deg');

ctrl.scaleThreeThrCtlr;

%% simulate
try
    set_param('mainModel','SimulationCommand','update');
catch
    set_param('mainModel','SimulationCommand','update');
end
  
simWithMonitor('mainModel',2);


%% post process
run_no = 1;
if lengthScale ~= 1 || densityScale ~=1
    run_no = 2;
end
scaledModel_postProcess


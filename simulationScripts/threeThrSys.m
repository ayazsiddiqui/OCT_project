% clear it all
clear all
clear mex;
fclose all;
clear

clc
format compact
% close all

%% simtime
plot_animation = 0;
make_video = 0;

%% common parameters
lengthScale = 1/1;
densityScale = 1/1;
numTethers = 3;
thrNumNodes = 2;
numTurbines = 2;

sim_time = 600*sqrt(lengthScale);

%% set variants
vhcl_variant = 'partitionedLiftingBodyVariant';
thr_variant = 'KelvinVoigtTetherVariant';
wnch_variant = 'PureIntegratorWinchVariant';
gnd_variant = 'FixedGrounStationVariant';

dts = 0.05*sqrt(lengthScale);
tVec = 0:dts:sim_time;

%% setpoints

% altitude
altitudeSP = 50*ones(size(tVec)).*lengthScale;

% pitch
pitchSP = 7*(pi/180)*ones(size(tVec));

% roll
Yswitch = 10*lengthScale;
rollAmp = 20;
rollPeriod = 100*sqrt(lengthScale);
startRoll = 0;
rollSP = (pi/180)*rollAmp*sign(sin((2*pi/rollPeriod)*tVec));

% yaw
yawSP = 0*(pi/180)*ones(size(tVec));

altitudeSP = timeseries(altitudeSP,tVec);
pitchSP = timeseries(pitchSP,tVec);
rollSP = timeseries(rollSP,tVec);
yawSP = timeseries(yawSP,tVec);


%% environment
env = ENV.environment;

env.setFluidDensity(1000,'kg/m^3');
env.setInertialFlowVel([1 0 0],'m/s');

%% lifiting body
vhcl = PLT.vehicle;

vhcl.setNumTethers(numTethers,'');
vhcl.setNumTurbines(numTurbines,'');
vhcl.setBuoyFactor(0.9,'');
vhcl.setFluidDensity(env.fluidDensity.Value,'kg/m^3');

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
vhcl.setRwingLE_cm([-1;0;0],'m');
vhcl.setWingChord(1,'m');
vhcl.setWingAR(10,'');
vhcl.setWingTR(0.8,'');
vhcl.setWingSweep(2,'deg');
vhcl.setWingDihedral(0,'deg');
vhcl.setWingIncidence(0,'deg');
vhcl.setWingNACA('4412','');
vhcl.setWingClMax(1.75,'');
vhcl.setWingClMin(-1.75,'');

% % % H-stab
vhcl.setRhsLE_wingLE([6;0;0],'m');
vhcl.setHsChord(0.6,'m');
vhcl.setHsAR(8,'');
vhcl.setHsTR(0.8,'');
vhcl.setHsSweep(5,'deg');
vhcl.setHsDihedral(0,'deg');
vhcl.setHsIncidence(0,'deg');
vhcl.setHsNACA('0012','');
vhcl.setHsClMaxl(1.75,'');
vhcl.setHsClMin(-1.75,'');

% % % V-stab
vhcl.setRvs_wingLE([6;0;0],'m');
vhcl.setVsChord(0.6,'m');
vhcl.setVsSpan(2.5,'m');
vhcl.setVsTR(0.8,'');
vhcl.setVsSweep(10,'deg');
vhcl.setVsNACA('0012','');
vhcl.setVsClMax(1.75,'');
vhcl.setVsClMin(-1.75,'');

% % % initial conditions
vhcl.setInitialCmPos([0;0;50],'m');
vhcl.setInitialCmVel([0;0;0],'m/s');
vhcl.setInitialEuler([0;1;0]*pi/180,'rad');
vhcl.setInitialAngVel([0;0;0],'rad/s');

% % % data file name
vhcl.setFluidCoeffsFileName('someFile','');

% % % load/generate fluid dynamic data
vhcl.calcFluidDynamicCoefffs

% % % calc aerodynamic center
% vhcl.calcDesginFluidDynamicCenter

% % % plot
% vhcl.plot
% vhcl.plotCoeffPolars

%% turbines
turb = PLT.turbine;

turb.setNumTurbines(numTurbines,'');
turb.setTurbDiameter(0*ones(1,numTurbines),'m')
turb.setTurbDragCoeff(0.8*ones(1,numTurbines),'');
turb.setTurbPowerCoeff(0.5*ones(1,numTurbines),'');

%% ground station
gnd = PLT.gndStn;

gnd.setNumTethers(numTethers,'');

gnd.setIzz(100,'kg*m^2');
gnd.setDampingCoeff(10,'N*m*s');
gnd.setFreeSpinSwitch(0,'');

gnd.setThrAttchPts(vhcl.thrAttchPts.Value.*[ones(2,numTethers);zeros(1,numTethers)],'m');

% % % initial conditions
gnd.setInitialEuler(0,'rad');
gnd.setInitialAngVel(0,'rad/s');


%% tethers
thr = PLT.tether;

thr.setNumTethers(numTethers,'');

thr.setNumNodes(thrNumNodes,'');
thr.setThrDensity(1300*ones(1,numTethers),'kg/m^3');
thr.setThrYoungs(3.8e9*ones(1,numTethers),'N/m^2');
thr.setThrDampingRatio(0.05*ones(1,numTethers),'');
thr.setThrDragCoeff(0.5*ones(1,numTethers),'');

% design tether
maxAppFlowMultiplier = 2;
maxPercentageElongation = 0.05;
valThrDia = thr.recommendTetherDiameter(vhcl,env,maxAppFlowMultiplier,maxPercentageElongation);

% set tether diameter
thr.setThrDiameter(valThrDia,'m');

%% winches
wnch = PLT.winch;

wnch.setNumTethers(numTethers,'');

wnch.setWnchMaxTugSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchMaxReleaseSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchTimeConstant(0.05*ones(1,numTethers),'s');

valInitLength = wnch.recommendInitTetherLength(vhcl,gnd,thr,env);
wnch.setInitThrLength(valInitLength,'m');


%% controller
ctrl = CTR.threeThrCtlr;

ctrl.setNumTethers(numTethers,'');
% altitude tether control gains
ctrl.setAltiTetherKp(0.0,'(m/s)/(m)')
ctrl.setAltiTetherKi(0,'(m/s)/(m*s)')
ctrl.setAltiTetherKd(0,'(m/s)/(m/s)')
ctrl.setAltiTetherTau(1,'s')
ctrl.setAltiErrorSat(5,'m')

% pitch tether control gains
ctrl.setPitchTetherKp(1*2,'(m/s)/(rad)')
ctrl.setPitchTetherKi(0,'(m/s)/(rad*s)')
ctrl.setPitchTetherKd(1*4,'(m/s)/(rad/s)')
ctrl.setPitchTetherTau(0.1,'s')

% roll tether control gains
ctrl.setRollTetherKp(1*4,'(m/s)/(rad)')
ctrl.setRollTetherKi(0,'(m/s)/(rad*s)')
ctrl.setRollTetherKd(1*12,'(m/s)/(rad/s)')
ctrl.setRollTetherTau(0.01,'s')

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


%% scale
env.scale(lengthScale,densityScale);
vhcl.scale(lengthScale,densityScale);
turb.scale(lengthScale,densityScale);
gnd.scale(lengthScale,densityScale);
thr.scale(lengthScale,densityScale);
wnch.scale(lengthScale,densityScale);
ctrl.scale(lengthScale,densityScale);

%% simulate
try
    %     open_system('mainModel');
    simWithMonitor('mainModel',2);
catch
    %     open_system('mainModel');
    simWithMonitor('mainModel',2);
end


%% post process
% scaledModel_postProcess

if plot_animation == 1
    fullKitePlotter
end


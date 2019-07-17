clear
clc
format compact

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

%% controller
ctrl = CTR.oneThrCtlr;

ctrl.setLengthScale(lengthScale,'');
ctrl.setDensityScale(densityScale,'');
ctrl.setNumTethers(numTethers,'');

% aileron gains
ctrl.setAileronKp(1,'deg/deg');
ctrl.setAileronKi(0,'deg/(deg*s)');
ctrl.setAileronKd(0,'deg/(deg/s)');
ctrl.setAileronTau(0.5,'s');
ctrl.setAileronMaxDef(30,'deg');

% elevator gains
ctrl.setElevatorKp(1,'deg/deg');
ctrl.setElevatorKi(0,'deg/(deg*s)');
ctrl.setElevatorKd(0,'deg/(deg/s)');
ctrl.setElevatorTau(0.5,'s');
ctrl.setElevatorMaxDef(30,'deg');

% rudder gains
ctrl.setRudderKp(1,'deg/deg');
ctrl.setRudderKi(0,'deg/(deg*s)');
ctrl.setRudderKd(0,'deg/(deg/s)');
ctrl.setRudderTau(0.5,'s');
ctrl.setRudderMaxDef(30,'deg');

ctrl.scaleThreeThrCtlr;


%% test signal
alti_sp = 0;
pitch_sp = 0*pi/180;
roll_sp = 0*pi/180;
yaw_sp = 5*pi/180;

Rcm_o = [0;0;0];
Vcm_o = zeros(3,1);
eulerAng = [0;0;0]*pi/180;
OwB = [0;0;0];

flow = [0;0;0];

open_system('oneThrCtrl_th')
sim('oneThrCtrl_th')










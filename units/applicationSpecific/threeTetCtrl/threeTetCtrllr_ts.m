clear
clc
format compact

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

%% controller
ctrl = CTR.threeThrCtlr;

ctrl.setLengthScale(lengthScale,'');
ctrl.setDensityScale(densityScale,'');
ctrl.setNumTethers(numTethers,'');

% altitude tether control gains
ctrl.setAltiTetherKp(1,'(m/s)/m')
ctrl.setAltiTetherKi(0,'(m/s)/(m*s)')
ctrl.setAltiTetherKd(0,'(m/s)/(m/s)')
ctrl.setAltiTetherTau(1,'s')
ctrl.setAltiErrorSat(5,'m')

% pitch tether control gains
ctrl.setPitchTetherKp(1,'(m/s)/rad')
ctrl.setPitchTetherKi(0,'(m/s)/(rad*s)')
ctrl.setPitchTetherKd(0,'(m/s)/(rad/s)')
ctrl.setPitchTetherTau(1,'s')

% roll tether control gains
ctrl.setRollTetherKp(1,'(m/s)/rad')
ctrl.setRollTetherKi(0,'(m/s)/(rad*s)')
ctrl.setRollTetherKd(0,'(m/s)/(rad/s)')
ctrl.setRollTetherTau(1,'s')

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
pitch_sp = 1*pi/180;
roll_sp = 0*pi/180;
yaw_sp = 0*pi/180;

Rcm_o = [0;0;0];
Vcm_o = zeros(3,1);
eulerAng = [0;0;0]*pi/180;
OwB = [0;0;0];

flow = [0;0;0];

open_system('threeTetCtrllr_th')
sim('threeTetCtrllr_th')










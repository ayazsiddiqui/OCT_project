clear
clc
format compact

% this is the build script for creating a controller using class definition
% 'threeThrCtlr' for a three tethered system that is being used by ayaz

% the script saves the variables 'ctrl' and to a .mat file

%% controller
numTethers = 3;

ctrl = CTR.threeThrCtlr;

ctrl.setNumTethers(numTethers,'');
% altitude tether control gains
ctrl.setAltiTetherKp(0.0,'(m/s)/(m)')
ctrl.setAltiTetherKi(0,'(m/s)/(m*s)')
ctrl.setAltiTetherKd(0,'(m/s)/(m/s)')
ctrl.setAltiTetherTau(1,'s')
ctrl.setAltiErrorSat(5,'m')

% pitch tether control gains
ctrl.setPitchTetherKp(1*0,'(m/s)/(rad)')
ctrl.setPitchTetherKi(0,'(m/s)/(rad*s)')
ctrl.setPitchTetherKd(1*0,'(m/s)/(rad/s)')
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

%% save file in its respective directory
saveBuildFile(ctrl,'ctrl',mfilename);



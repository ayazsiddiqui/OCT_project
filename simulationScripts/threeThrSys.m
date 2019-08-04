% clear it all
clear 
clc
format compact
close all

cd(fileparts(mfilename('fullpath')));

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
load('threeThrEnv.mat');

env.setInertialFlowVel([1 0 0],'m/s');

%% vehicle
load('threeThrVhcl.mat');

% % % initial conditions
vhcl.setInitialCmPos([0;0;50],'m');
vhcl.setInitialCmVel([0;0;0],'m/s');
vhcl.setInitialEuler([0;1;0]*pi/180,'rad');
vhcl.setInitialAngVel([0;0;0],'rad/s');

% % % calc aerodynamic center
% vhcl.calcDesginFluidDynamicCenter

% % % plot
% vhcl.plot
% vhcl.plotCoeffPolars

%% turbine
load('dummyTurbines.mat');

%% ground station
load('threeThrGnd.mat')

% % % initial conditions
gnd.setInitialEuler(0,'rad');
gnd.setInitialAngVel(0,'rad/s');

%% tethers
load('threeThrtethers.mat');

%% winches
load('threeThrWnch.mat');

valInitLength = wnch.recommendInitTetherLength(vhcl,gnd,thr,env);
wnch.setInitThrLength(valInitLength,'m');


%% controller
load('threeThrCtrl.mat');

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
        open_system('mainModel');
    simWithMonitor('mainModel',2);
catch
    %     open_system('mainModel');
    simWithMonitor('mainModel',2);
end


%% post process
scaledModel_postProcess

if plot_animation == 1
    fullKitePlotter
end


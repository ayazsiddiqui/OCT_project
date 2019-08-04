clear
clc
format compact

% this is the build script for creating a environment using class definition
% 'environment' for a three tethered system that is being used by ayaz

% the script saves the variables 'env' and to a .mat file

%% environment
env = ENV.environment;

env.setFluidDensity(1000,'kg/m^3');
env.setInertialFlowVel([1 0 0],'m/s');

%% save file in its respective directory
saveBuildFile(env,'env',mfilename);

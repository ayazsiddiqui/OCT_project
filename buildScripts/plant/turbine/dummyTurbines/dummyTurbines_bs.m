clear
clc
format compact

% this is the build script for creating a turbines using class definition
% 'turbine' for a three tethered system that is being used by ayaz

% the script saves the variables 'turb' and to a .mat file

%% turbines
numTurbines = 2;


turb = PLT.turbine;

turb.setNumTurbines(numTurbines,'');
turb.setTurbDiameter(0*ones(1,numTurbines),'m')
turb.setTurbDragCoeff(0.8*ones(1,numTurbines),'');
turb.setTurbPowerCoeff(0.5*ones(1,numTurbines),'');

%% save file in its respective directory
saveBuildFile(turb,'turb',mfilename);



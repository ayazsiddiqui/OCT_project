clear
clc
format compact

% this is the build script for creating a tether using class definition
% 'tether' for a three tethered system that is being used by ayaz

% the script saves the variables 'thr' and to a .mat file

%% tether variant
thr_variant = 'KelvinVoigtTetherVariant';

%% tethers
numTethers = 3;
thrNumNodes = 2;

thr = PLT.tether;

thr.setNumTethers(numTethers,'');

thr.setNumNodes(thrNumNodes,'');
thr.setThrDensity(1300*ones(1,numTethers),'kg/m^3');
thr.setThrYoungs(4e9*ones(1,numTethers),'N/m^2');
thr.setThrDampingRatio(0.05*ones(1,numTethers),'');
thr.setThrDragCoeff(0.5*ones(1,numTethers),'');

% design tether
maxAppFlowMultiplier = 2;
maxPercentageElongation = 0.05;

thrDia = 0.0075;

% set tether diameter
thr.setThrDiameter(thrDia.*[1 sqrt(2) 1],'m');

%% save file in its respective directory
saveBuildFile(thr,'thr',mfilename,'variantVariableName','thr_variant','variantVariable',thr_variant);


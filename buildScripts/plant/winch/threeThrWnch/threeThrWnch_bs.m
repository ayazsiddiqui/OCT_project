clear
clc
format compact

% this is the build script for creating a winches using class definition
% 'winch' for a three tethered system that is being used by ayaz

% the script saves the variables 'wnch' and to a .mat file

%% winch variant
wnch_variant = 'PureIntegratorWinchVariant';

%% winches
numTethers = 3;

wnch = PLT.winch;

wnch.setNumTethers(numTethers,'');

wnch.setWnchMaxTugSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchMaxReleaseSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchTimeConstant(0.05*ones(1,numTethers),'s');

%% save file in its respective directory
saveBuildFile(wnch,'wnch',mfilename,'variantVariableName','wnch_variant','variantVariable',wnch_variant);



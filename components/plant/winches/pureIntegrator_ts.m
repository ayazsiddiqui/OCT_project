clear
clc
format compact

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

%% winches
wnch = PLT.winch;

wnch.setLengthScale(lengthScale,'');
wnch.setDensityScale(densityScale,'');
wnch.setNumTethers(numTethers,'');

wnch.setWnchMaxTugSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchMaxReleaseSpeed(1*ones(1,numTethers),'m/s');
wnch.setWnchTimeConstant(1*ones(1,numTethers),'s');

wnch.setInitThrLength(100*ones(1,numTethers),'m');

wnch.scaleWinch;

%% test signals
testRelease = -4*ones(1,numTethers);




clear
clc
format compact

% this is the build script for creating a ground station using class definition
% 'gndStn' for a three tethered system that is being used by ayaz

% the script saves the variables 'gnd' and 'gnd_variant' to a .mat file

%% variant
gnd_variant = 'FixedGrounStationVariant';

%% ground station
numTethers = 3;

gnd = PLT.gndStn;

gnd.setNumTethers(numTethers,'');

gnd.setIzz(100,'kg*m^2');
gnd.setDampingCoeff(10,'N*m*s');
gnd.setFreeSpinSwitch(0,'');

gnd.setThrAttchPts([0.8254   -5.0000         0;
    5.6000  0 0;
    -0.8254    5.0000 0]'.*[ones(2,numTethers);zeros(1,numTethers)],'m');


%% save file in its respective directory
saveBuildFile(gnd,'gnd',mfilename,'variantVariableName','gnd_variant','variantVariable',gnd_variant);



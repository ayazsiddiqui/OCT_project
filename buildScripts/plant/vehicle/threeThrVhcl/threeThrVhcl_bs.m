clear
clc
format compact

% this is the build script for creating a vechile using class definition
% 'vehicle' for a three tethered system that is being used by ayaz

% the script saves the variables 'vhcl' and 'vhcl_variant' to a .mat file

%% set variant
vhcl_variant = 'partitionedLiftingBodyVariant';

%% lifiting body
numTethers = 3;
numTurbines = 2;

vhcl = PLT.vehicle;

vhcl.setNumTethers(numTethers,'');
vhcl.setNumTurbines(numTurbines,'');
vhcl.setBuoyFactor(0.9,'');
vhcl.setFluidDensity(1000,'kg/m^3');

% % % volume and inertias
vhcl.setVolume(945352023.474*1e-9,'m^3');
vhcl.setIxx(6.303080401918E+09*1e-6,'kg*m^2');
vhcl.setIyy(2080666338.077*1e-6,'kg*m^2');
vhcl.setIzz(8.320369733598E+09*1e-6,'kg*m^2');
vhcl.setIxy(0,'kg*m^2');
vhcl.setIxz(81875397.942*1e-6,'kg*m^2');
vhcl.setIyz(0,'kg*m^2');
vhcl.setRcb_cm([0;0;0],'m');

% % % wing
vhcl.setRwingLE_cm([-1;0;0],'m');
vhcl.setWingChord(1,'m');
vhcl.setWingAR(10,'');
vhcl.setWingTR(0.8,'');
vhcl.setWingSweep(2,'deg');
vhcl.setWingDihedral(0,'deg');
vhcl.setWingIncidence(0,'deg');
vhcl.setWingNACA('4412','');
vhcl.setWingClMax(1.75,'');
vhcl.setWingClMin(-1.75,'');

% % % H-stab
vhcl.setRhsLE_wingLE([6;0;0],'m');
vhcl.setHsChord(0.6,'m');
vhcl.setHsAR(8,'');
vhcl.setHsTR(0.8,'');
vhcl.setHsSweep(5,'deg');
vhcl.setHsDihedral(0,'deg');
vhcl.setHsIncidence(0,'deg');
vhcl.setHsNACA('0012','');
vhcl.setHsClMaxl(1.75,'');
vhcl.setHsClMin(-1.75,'');

% % % V-stab
vhcl.setRvs_wingLE([6;0;0],'m');
vhcl.setVsChord(0.6,'m');
vhcl.setVsSpan(2.5,'m');
vhcl.setVsTR(0.8,'');
vhcl.setVsSweep(10,'deg');
vhcl.setVsNACA('0012','');
vhcl.setVsClMax(1.75,'');
vhcl.setVsClMin(-1.75,'');

% % % data file name
vhcl.setFluidCoeffsFileName('someFile','');

% % % load/generate fluid dynamic data
vhcl.calcFluidDynamicCoefffs

%% save file in its respective directory
saveBuildFile(vhcl,'vhcl',mfilename,'variantVariableName','vhcl_variant','variantVariable',vhcl_variant);



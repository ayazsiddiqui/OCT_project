clear
clc
format compact

% create class instance
tp = plant(3,2);

tp.ScaleFactor = 1;

% set vehicle values
tp.vehicle.mass.value = 8.9360e+04;
tp.vehicle.MI.value = [14330000 0 0; 0 143200 0; 0 0 15300000];
tp.vehicle.volume.value = 111.7000;
tp.vehicle.Rcb_cm.value = [0;0;0];
tp.vehicle.Rcm_wingLE.value = [0;0;0];

% attachment points
tp.vehicle.tetherAttchPts(1).value = [-2.5000; -20.0000; -0.3750];
tp.vehicle.tetherAttchPts(2).value = [21.2500; 0; -0.3750];
tp.vehicle.tetherAttchPts(3).value = [-2.5000; 20.0000; -0.3750];

% turbines
tp.turbines(1).Rturb_cm.value = [2.5000; -20.2500; 0];
tp.turbines(1).diameter.value = 8.7;
tp.turbines(1).powerCoeff.value = 0.5;
tp.turbines(1).dragCoeff.value = 0.8;

tp.turbines(2).Rturb_cm.value = [2.5000; 20.2500; 0];
tp.turbines(2).diameter.value = 8.7;
tp.turbines(2).powerCoeff.value = 0.5;
tp.turbines(2).dragCoeff.value = 0.8;

% tethers
tp.tethers(1).numNodes.value = 4;
tp.tethers(1).diameter.value = 0.055;
tp.tethers(1).youngsModulus.value = 3.8e9;
tp.tethers(1).dampingRation.value = 0.05;
tp.tethers(1).dragCoeff.value = 0.5;
tp.tethers(1).density.value = 1300;






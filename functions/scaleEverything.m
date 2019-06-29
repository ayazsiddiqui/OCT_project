function [plant,environment,controller,simTime] = scaleEverything(...
    plant,environment,controller,simTime,altitudeSetpoint)

%% dummy variables
% length scale factor
Lscale = plant.lengthScaleFactor;
% density scale factor
Dscale = plant.densityScaleFactor;

% others
pl = plant;
en = environment;
ct = controller;

%% scale simulation time
simTime = simTime*sqrt(Lscale);

%% scale setpoint
altitudeSetpoint.Data = altitudeSetpoint.Data*Lscale;
altitudeSetpoint.Time = altitudeSetpoint.Time*sqrt(Lscale);


%% scale environment
en.flowDensity.value = en.flowDensity.value*(Dscale);
en.iniertialFlowVel.value = en.iniertialFlowVel.value*sqrt(Lscale);


%% Scale plant

% vehicle
pl.vehicle.mass.value = pl.vehicle.mass.value*(Lscale^3);
pl.vehicle.MI.value = pl.vehicle.MI.value*(Lscale^5);
pl.vehicle.volume.value = pl.vehicle.volume.value*(Lscale^3);
pl.vehicle.Rcb_cm.value = pl.vehicle.Rcb_cm.value*(Lscale);
pl.vehicle.Rcm_wingLE.value = pl.vehicle.Rcm_wingLE.value*(Lscale);
pl.vehicle.ini_Rcm_o.value = pl.vehicle.ini_Rcm_o.value*(Lscale);
pl.vehicle.ini_O_Vcm_o.value = pl.vehicle.ini_O_Vcm_o.value*sqrt(Lscale);
pl.vehicle.ini_OwB.value = pl.vehicle.ini_OwB.value*(1/sqrt(Lscale));

pl = pl.calcAddedMass(en);

% turbine
for ii = 1:length(pl.turbines)
    pl.turbines(ii).Rturb_cm.value = pl.turbines(ii).Rturb_cm.value*(Lscale);
    pl.turbines(ii).diameter.value = pl.turbines(ii).diameter.value*(Lscale);
end

% tethers
for ii = 1:length(pl.tethers)
    pl.tethers(ii).diameter = pl.tethers(ii).diameter*(Lscale)*(Dscale^(1/1.985));
    pl.tethers(ii).youngsModulus = pl.tethers(ii).youngsModulus*(Lscale);
    pl.tethers(ii).density = pl.tethers(ii).density*(Dscale);
    
end

% ground station
pl.gndStation.Izz.value = pl.gndStation.Izz.value*(Lscale^5);
pl.gndStation.dampCoeff.value = pl.gndStation.dampCoeff.value*(Lscale^4.5);
pl.gndStation.ini_platform_vel.value = pl.gndStation.ini_platform_vel.value*(1/sqrt(Lscale));

% winches
for ii = 1:length(pl.winches)
    pl.winches(ii).maxSpeed = pl.winches(ii).maxSpeed*(sqrt(Lscale));
    pl.winches(ii).timeConstant = pl.winches(ii).timeConstant*(1/sqrt(Lscale));
    pl.winches(ii).initTetherLength = pl.winches(ii).initTetherLength*(Lscale);
end   
    
%% Scale controller

% tether controllers
ct.tethers.altiTetherKp = ct.tethers.altiTetherKp*(sqrt(Lscale));
ct.tethers.altiTetherKi = ct.tethers.altiTetherKi*(1/Lscale);
ct.tethers.altiTetherKd = ct.tethers.altiTetherKd*1;
ct.tethers.altiTetherTau = ct.tethers.altiTetherTau*(1/sqrt(Lscale));

ct.tethers.pitchTetherKp = ct.tethers.pitchTetherKp*(sqrt(Lscale));
ct.tethers.pitchTetherKi = ct.tethers.pitchTetherKi*1;
ct.tethers.pitchTetherKd = ct.tethers.pitchTetherKd*(Lscale);
ct.tethers.pitchTetherTau = ct.tethers.pitchTetherTau*(1/sqrt(Lscale));

ct.tethers.rollTetherKp = ct.tethers.rollTetherKp*(sqrt(Lscale));
ct.tethers.rollTetherKi = ct.tethers.rollTetherKi*1;
ct.tethers.rollTetherKd = ct.tethers.rollTetherKd*(Lscale);
ct.tethers.rollTetherTau = ct.tethers.rollTetherTau*(1/sqrt(Lscale));

% control surface controllers
ct.controlSurfaces.aileronKp = ct.controlSurfaces.aileronKp;
ct.controlSurfaces.aileronKi = ct.controlSurfaces.aileronKi*(1/sqrt(Lscale));
ct.controlSurfaces.aileronKd = ct.controlSurfaces.aileronKd*(sqrt(Lscale));
ct.controlSurfaces.aileronTau = ct.controlSurfaces.aileronTau*(1/sqrt(Lscale));

ct.controlSurfaces.elevatorKp = ct.controlSurfaces.elevatorKp;
ct.controlSurfaces.elevatorKi = ct.controlSurfaces.elevatorKi*(1/sqrt(Lscale));
ct.controlSurfaces.elevatorKd = ct.controlSurfaces.elevatorKd*(sqrt(Lscale));
ct.controlSurfaces.elevatorTau = ct.controlSurfaces.elevatorTau*(1/sqrt(Lscale));

% that should work

end
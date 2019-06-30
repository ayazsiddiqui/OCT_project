function [s_plant,s_environment,s_controller,s_simTime,altitudeSetpoint,pitchSP,rollSP]...
    = scaleEverything(scaleFactors,plant,environment,controller,simTime,altitudeSetpoint,pitchSP,rollSP)

%% dummy variables
% length scale factor
Lscale = scaleFactors(1);
% density scale factor
Dscale = scaleFactors(2);

% others
o_pl = plant;
o_en = environment;
o_ct = controller;

% create dummy classdefs
n_pl = sysParam.plant_v2(o_pl.numTethers,o_pl.numTurbines);
n_en = sysParam.env;

%% scale simulation time
s_simTime = simTime*(Lscale^0.5);

%% scale setpoint
altitudeSetpoint.Data = altitudeSetpoint.Data.*(Lscale);
altitudeSetpoint.Time = altitudeSetpoint.Time.*(Lscale^0.5);

pitchSP.Time = pitchSP.Time.*(Lscale^0.5);
rollSP.Time = rollSP.Time.*(Lscale^0.5);


%% scale environment
n_en.gravAccel.value = 9.81;
n_en.flowDensity.value = o_en.flowDensity.value*(Dscale);
n_en.inertialFlowVel.value = o_en.inertialFlowVel.value*(Lscale^0.5);

n_pl.lengthScaleFactor = Lscale;
n_pl.densityScaleFactor = Dscale;
n_pl.buoyancyFactor = o_pl.buoyancyFactor;


%% Scale plant

% vehicle
n_pl.vehicle.MI.value = o_pl.vehicle.MI.value.*(Lscale^5);
n_pl.vehicle.volume.value = o_pl.vehicle.volume.value*(Lscale^3);
n_pl.vehicle.mass.value = (1/n_pl.buoyancyFactor)*n_pl.vehicle.volume.value*n_en.flowDensity.value;
n_pl.vehicle.Rcb_cm.value = o_pl.vehicle.Rcb_cm.value*(Lscale);
n_pl.vehicle.Rcm_wingLE.value = o_pl.vehicle.Rcm_wingLE.value*(Lscale);
n_pl.vehicle.ini_Rcm_o.value = o_pl.vehicle.ini_Rcm_o.value*(Lscale);
n_pl.vehicle.ini_O_Vcm_o.value = o_pl.vehicle.ini_O_Vcm_o.value*(Lscale^0.5);
n_pl.vehicle.ini_euler.value = o_pl.vehicle.ini_euler.value;
n_pl.vehicle.ini_OwB.value = o_pl.vehicle.ini_OwB.value*(1/(Lscale^0.5));

n_pl.aeroDataFileName = o_pl.aeroDataFileName;

n_pl = n_pl.calcAddedMass(n_en);

% turbine
for ii = 1:length(n_pl.turbines)
    n_pl.turbines(ii).Rturb_cm.value = o_pl.turbines(ii).Rturb_cm.value*(Lscale);
    n_pl.turbines(ii).diameter.value = o_pl.turbines(ii).diameter.value*(Lscale);
    n_pl.turbines(ii).powerCoeff.value = o_pl.turbines(ii).powerCoeff.value;
    n_pl.turbines(ii).dragCoeff.value = o_pl.turbines(ii).dragCoeff.value;

end

% tethers
for ii = 1:length(n_pl.tethers)
    n_pl.tethers(ii).numNodes = o_pl.tethers(ii).numNodes;
    n_pl.tethers(ii).diameter = o_pl.tethers(ii).diameter*(Lscale)*(Dscale^(1/1.985));
    n_pl.tethers(ii).youngsModulus = o_pl.tethers(ii).youngsModulus*(Lscale);
    n_pl.tethers(ii).dampingRatio = o_pl.tethers(ii).dampingRatio;
    n_pl.tethers(ii).dragCoeff = o_pl.tethers(ii).dragCoeff;
    n_pl.tethers(ii).density = o_pl.tethers(ii).density*(Dscale);
    
end

% maxAppFlowMultiplier = 4;
% maxPercentageElongation = 0.01;
% n_pl = n_pl.designTetherDiameter(n_en,maxAppFlowMultiplier,maxPercentageElongation);

% ground station
n_pl.gndStation.rotationSwitch.value = o_pl.gndStation.rotationSwitch.value;

n_pl.gndStation.Izz.value = o_pl.gndStation.Izz.value*(Lscale^5);
n_pl.gndStation.dampCoeff.value = o_pl.gndStation.dampCoeff.value*(Lscale^4.5);

n_pl.gndStation.ini_platform_ang.value = o_pl.gndStation.ini_platform_ang.value;
n_pl.gndStation.ini_platform_vel.value = o_pl.gndStation.ini_platform_vel.value*(1/(Lscale^0.5));

% winches
for ii = 1:length(n_pl.winches)
    n_pl.winches(ii).maxSpeed = o_pl.winches(ii).maxSpeed*(Lscale^0.5);
    n_pl.winches(ii).timeConstant = o_pl.winches(ii).timeConstant*(Lscale^0.5);
    n_pl.winches(ii).initTetherLength = o_pl.winches(ii).initTetherLength*(Lscale);
end

% n_pl = n_pl.setTetherInitLength(n_en);

    
%% Scale controller
n_ct.tethers.transformMat = o_ct.tethers.transformMat;

% tether controllers
n_ct.tethers.altiTetherKp = o_ct.tethers.altiTetherKp*(1/Lscale^0.5);
n_ct.tethers.altiTetherKi = o_ct.tethers.altiTetherKi*(1/Lscale);
n_ct.tethers.altiTetherKd = o_ct.tethers.altiTetherKd;
n_ct.tethers.altiTetherTau = o_ct.tethers.altiTetherTau*(Lscale^0.5);
n_ct.tethers.altiErrorSat = o_ct.tethers.altiErrorSat*(1/Lscale);

n_ct.tethers.pitchTetherKp = o_ct.tethers.pitchTetherKp*(Lscale^0.5);
n_ct.tethers.pitchTetherKi = o_ct.tethers.pitchTetherKi;
n_ct.tethers.pitchTetherKd = o_ct.tethers.pitchTetherKd*(Lscale);
n_ct.tethers.pitchTetherTau = o_ct.tethers.pitchTetherTau*(Lscale^0.5);

n_ct.tethers.rollTetherKp = o_ct.tethers.rollTetherKp*(Lscale^0.5);
n_ct.tethers.rollTetherKi = o_ct.tethers.rollTetherKi;
n_ct.tethers.rollTetherKd = o_ct.tethers.rollTetherKd*(Lscale);
n_ct.tethers.rollTetherTau = o_ct.tethers.rollTetherTau*(Lscale^0.5);

% control surface controllers
n_ct.controlSurfaces.aileronKp = o_ct.controlSurfaces.aileronKp;
n_ct.controlSurfaces.aileronKi = o_ct.controlSurfaces.aileronKi*(1/(Lscale^0.5));
n_ct.controlSurfaces.aileronKd = o_ct.controlSurfaces.aileronKd*(Lscale^0.5);
n_ct.controlSurfaces.aileronTau = o_ct.controlSurfaces.aileronTau*(Lscale^0.5);
n_ct.controlSurfaces.aileronMaxDef = o_ct.controlSurfaces.aileronMaxDef;


n_ct.controlSurfaces.elevatorKp = o_ct.controlSurfaces.elevatorKp;
n_ct.controlSurfaces.elevatorKi = o_ct.controlSurfaces.elevatorKi*(1/(Lscale^0.5));
n_ct.controlSurfaces.elevatorKd = o_ct.controlSurfaces.elevatorKd*(Lscale^0.5);
n_ct.controlSurfaces.elevatorTau = o_ct.controlSurfaces.elevatorTau*(Lscale^0.5);
n_ct.controlSurfaces.elevatorMaxDef = o_ct.controlSurfaces.elevatorMaxDef;


% that should work
s_environment = n_en;
s_plant = n_pl;
s_controller = n_ct;




end
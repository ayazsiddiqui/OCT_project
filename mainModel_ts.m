clear
clc
format compact

%% environment conditions
flowVel = [0.2;0;0];
rho_fluid = 1000;
grav = 9.81;

%% lifting body parameters
% mass and buoyancy parameters
mass = 8.9360e+04;
m_added = [1.8017e+04 0 0;0 1.5825e+06 0;0 0 8.0374e+05];
MI = [14330000 0 0; 0 143200 0; 0 0 15300000];
vol = 111.7000;
Rcb_cm = [0;0;0];

class_liftingBdy(1).tetherAttchPt = [-2.5000; -20.0000; -0.3750];
class_liftingBdy(2).tetherAttchPt = [21.2500; 0; -0.3750];
class_liftingBdy(3).tetherAttchPt = [-2.5000; 20.0000; -0.3750];

class_liftingBdy = reshape(class_liftingBdy,1,[]);

% setpoints
altitudeSP = 200;
pitchSP = 2*pi/180;
rollSP = 0*pi/180;

% initial operating conditions
ini_Rcm_o = [0; 0; 200];
ini_O_Vcm_o = [0; 0; 0];
ini_euler = [0; 0.1222; 0];
ini_OwB = [0; 0; 0];

%% turbine parameters
% turbine struct
turbine(1).Rturb_cm = [2.5000; -20.2500; 0];
turbine(1).powerCoeff = 0.5;
turbine(1).dragCoeff = 0.8;
turbine(1).diameter =  8.7;

turbine(2).Rturb_cm = [2.5000; 20.2500; 0];
turbine(2).powerCoeff = 0.5;
turbine(2).dragCoeff = 0.8;
turbine(2).diameter =  8.7;

%% tether parameters
class_thr(1).numNodes = 4;
class_thr(1).diameter = 0.055;
class_thr(1).youngsModulus = 3.8e9;
class_thr(1).dampingRatio = 0.05;
class_thr(1).dragCoeff = 0.5;
class_thr(1).density = 1300;
class_thr(1).vehicleMass = 8.9360e+04;
class_thr(1).ini_R1_o = [-2.5000; -20.0000; 0.3750];
class_thr(1).ini_Rn_o = [-2.5271; -20.0000; 199.9325];

class_thr(2).numNodes = 4;
class_thr(2).diameter = 0.076;
class_thr(2).youngsModulus = 3.8e9;
class_thr(2).dampingRatio = 0.05;
class_thr(2).dragCoeff = 0.5;
class_thr(2).density = 1300;
class_thr(2).vehicleMass = 8.9360e+04;
class_thr(2).ini_R1_o = [21.2500; 0; 0.3750];
class_thr(2).ini_Rn_o = [21.0459; 0; 197.0381];

class_thr(3).numNodes = 4;
class_thr(3).diameter = 0.055;
class_thr(3).youngsModulus = 3.8e9;
class_thr(3).dampingRatio = 0.05;
class_thr(3).dragCoeff = 0.5;
class_thr(3).density = 1300;
class_thr(3).vehicleMass = 8.9360e+04;
class_thr(3).ini_R1_o = [-2.5000; 20.0000; 0.3750];
class_thr(3).ini_Rn_o = [-2.5271; 20.0000; 199.9325];

for ii = 1:length(class_thr)
    
    class_thr(ii).initNodePoss = [...
        linspace(class_thr(ii).ini_R1_o(1),class_thr(ii).ini_Rn_o(1),...
        class_thr(ii).numNodes);...
        linspace(class_thr(ii).ini_R1_o(2),class_thr(ii).ini_Rn_o(2),...
        class_thr(ii).numNodes);...
        linspace(class_thr(ii).ini_R1_o(3),class_thr(ii).ini_Rn_o(3),...
        class_thr(ii).numNodes)];
    
    class_thr(ii).initNodePoss = class_thr(ii).initNodePoss(:,2:end-1);
    class_thr(ii).initNodeVels = zeros(size(class_thr(ii).initNodePoss));
    
    n_ini_pos(:,:,ii) = class_thr(ii).initNodePoss;
    n_ini_vel(:,:,ii) = class_thr(ii).initNodeVels;
    
end

%% winch parameters
winchMaxSpeed = 0.4*ones(1,length(class_thr));
winchTimeConstant = ones(1,length(class_thr));
initTetherLength = [199.0076 196.0955 199.0076];

%% ground station parameters
Izz = 100;
c_damp = 10;

ini_platform_ang = 0;
ini_platform_vel = 0;

class_gndStn(1).tetherAttchPt = [-2.5000; -20.0000; 0.3750];
class_gndStn(2).tetherAttchPt = [21.2500; 0; 0.3750];
class_gndStn(3).tetherAttchPt = [-2.5000; 20.0000; 0.3750];

class_gndStn = reshape(class_gndStn,1,[]);

%% avl parameters
load('dsgnTest_1_lookupTables');

Sref = dsgnData.Sref;
Bref = dsgnData.Bref;
Cref = dsgnData.Cref;

%% define controller parameter
% tether command gains
ctrllr.tethers.transformMat = [1 .5 -.5; 1 -.5 0; 1 .5 .5];

ctrllr.tethers.altiTetherKp = 0.05;    % m/s per m
ctrllr.tethers.altiTetherKi = 0;
ctrllr.tethers.altiTetherKd = 5*ctrllr.tethers.altiTetherKp;
ctrllr.tethers.altiTetherTau = 31.8310;

ctrllr.tethers.pitchTetherKp = 1.5*0.5;   % m/s per rad
ctrllr.tethers.pitchTetherKi = 0;
ctrllr.tethers.pitchTetherKd = 2.5*ctrllr.tethers.pitchTetherKp;
ctrllr.tethers.pitchTetherTau = 0.7958;

ctrllr.tethers.rollTetherKp = 1*ctrllr.tethers.pitchTetherKp/(0.5*dsgnData.wing_AR);    % m/s per rad
ctrllr.tethers.rollTetherKi = 0;
ctrllr.tethers.rollTetherKd = 2*ctrllr.tethers.rollTetherKp;
ctrllr.tethers.rollTetherTau = 0.7958;

% control surface gains
ctrllr.controlSurfaces.aileronKp = 2;   % deg per deg
ctrllr.controlSurfaces.aileronKi = 0;
ctrllr.controlSurfaces.aileronKd = 2*ctrllr.controlSurfaces.aileronKp;
ctrllr.controlSurfaces.aileronTau = 0.2;
ctrllr.controlSurfaces.aileronMaxDef = 30;

ctrllr.controlSurfaces.elevatorKp = 1;  % deg per deg
ctrllr.controlSurfaces.elevatorKi = 0;
ctrllr.controlSurfaces.elevatorKd = 3*ctrllr.controlSurfaces.elevatorKp;
ctrllr.controlSurfaces.elevatorTau = 0.05;

ctrllr.controlSurfaces.elevatorMaxDef = 30;




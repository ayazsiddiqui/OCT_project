clear
clc
format compact

%% lifting body parameters
% mass and buoyancy parameters
mass = 8.9360e+04;
m_added = [1.8017e+04 0 0;0 1.5825e+06 0;0 0 8.0374e+05];
MI = [14330000 0 0; 0 143200 0; 0 0 15300000];
grav = 9.81;
vol = 111.7000;
rho_fluid = 1000;
Rcb_cm = [0;0;0];

% initial operating conditions
ini_Rcm_o = [0; 0; 200];
ini_O_Vcm_o = [0; 0; 0];
ini_euler = [0; 0.1222; 0];
ini_OwB = [0; 0; 0];

%% turbine parameters
% turbine struct
turbine(1).Rturb_cm = [2.5000; -20.2500; 0];
turbine(1).Cp = 0.5;
turbine(1).Cd = 0.8;
turbine(1).dia = 8.7;

turbine(2).Rturb_cm = [2.5000; 20.2500; 0];
turbine(2).Cp = 0.5;
turbine(2).Cd = 0.8;
turbine(2).dia = 8.7;

%% tether parameters
class_thr(1).R1_g = [-2.5000; -20.0000; 0.3750];
class_thr(1).Rn_cm = [-2.5000; -20.0000; -0.3750];
class_thr(1).numNodes = 4;
class_thr(1).diameter = 0.055;
class_thr(1).youngsModulus = 3.8e9;
class_thr(1).dampingRatio = 0.05;
class_thr(1).dragCoeff = 0.5;
class_thr(1).density = 1300;
class_thr(1).vehicleMass = 8.9360e+04;
class_thr(1).ini_R1_o = [-2.5000; -20.0000; 0.3750];
class_thr(1).ini_Rn_o = [-2.5271; -20.0000; 199.9325];

class_thr(2).R1_g = [21.2500; 0; 0.3750];
class_thr(2).Rn_cm = [21.2500; 0; -0.3750];
class_thr(2).numNodes = 4;
class_thr(2).diameter = 0.076;
class_thr(2).youngsModulus = 3.8e9;
class_thr(2).dampingRatio = 0.05;
class_thr(2).dragCoeff = 0.5;
class_thr(2).density = 1300;
class_thr(2).vehicleMass = 8.9360e+04;
class_thr(2).ini_R1_o = [21.2500; 0; 0.3750];
class_thr(2).ini_Rn_o = [21.0459; 0; 197.0381];

class_thr(3).R1_g = [-2.5000; 20.0000; 0.3750];
class_thr(3).Rn_cm = [-2.5000; 20.0000; -0.3750];
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
        linspace(class_thr(ii).ini_R1_o(1),class_thr(ii).ini_Rn_o(1),class_thr(ii).numNodes);...
        linspace(class_thr(ii).ini_R1_o(2),class_thr(ii).ini_Rn_o(2),class_thr(ii).numNodes);...
        linspace(class_thr(ii).ini_R1_o(3),class_thr(ii).ini_Rn_o(3),class_thr(ii).numNodes)];
    
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


%% avl parameters
load('dsgnTest_1_lookupTables');

Sref = dsgnData.Sref;
Bref = dsgnData.Bref;
Cref = dsgnData.Cref;

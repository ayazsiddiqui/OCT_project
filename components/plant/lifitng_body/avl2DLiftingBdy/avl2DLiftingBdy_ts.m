% clear
clear
clc
format compact

% start
%% mass and buoyancy parameters
mass = 10;
MI = eye(3);
grav = 9.81;
vol = 20;
rho_fluid = 1;
Rcb_cm = [0;0;0.1];

%% turbine parameters
% turbine struct
turbine(1).Rturb_cm = [0;10;0];
turbine(1).Cp = 0.5;
turbine(1).Cd = 0.8;
turbine(1).dia = 1;

turbine(2).Rturb_cm = [0;10;0];
turbine(2).Cp = 0.5;
turbine(2).Cd = 0.8;
turbine(2).dia = 1;

%% tether parameters
class_thr(1).R1_g = [0; -1; 0];
class_thr(1).Rn_cm = [0; -1; -1];

class_thr(2).R1_g = [1; 0; 0];
class_thr(2).Rn_cm = [1; 0; -1];

class_thr(3).R1_g = [0; 1; 0];
class_thr(3).Rn_cm = [0; 1; -1];

%% avl parameters
load('dsgnTest_1_lookupTables');

Sref = dsgnData.Sref;
Bref = dsgnData.Bref;
Cref = dsgnData.Cref;


%% testing signals
ini_Rcm_o = [0;0;0];
ini_O_Vcm_o = [0;0;0];
ini_euler = [0;0;0];
ini_OwB = [0;0;0];

euler = [0;0;0]*pi/180;
vFlow = [1;0;0];
vBdy = [0;0;0];







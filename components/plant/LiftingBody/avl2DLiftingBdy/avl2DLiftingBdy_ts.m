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
turbine(1).powerCoeff = 0.5;
turbine(1).dragCoeff = 0.8;
turbine(1).diameter = 1;

turbine(2).Rturb_cm = [0;10;0];
turbine(2).powerCoeff = 0.5;
turbine(2).dragCoeff = 0.8;
turbine(2).diameter = 1;

turbine = reshape(turbine,1,[]);

%% tether attachment points parameters
class_liftingBdy(1).tetherAttchPt = [0; -1; 0];
class_liftingBdy(2).tetherAttchPt = [1; 0; 0];
class_liftingBdy(3).tetherAttchPt = [0; 1; 0];

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







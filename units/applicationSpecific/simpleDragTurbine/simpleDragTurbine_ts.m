clear
clc
format compact

%% start
% turbine struct
turbine(1).Rturb_cm = [0;10;0];
turbine(1).Cp = 0.5;
turbine(1).Cd = 0.8;
turbine(1).dia = 1;

turbine(2).Rturb_cm = [0;10;0];
turbine(2).Cp = 0.5;
turbine(2).Cd = 0.8;
turbine(2).dia = 1;

% other params
V_flow = [1;0;0];
Vbody = [0;0;0];
euler = [0;0;0];

rho_fluid = 1;

sim('simpleDragTurbine_th')






clear
clc
format compact

%% define controller parameter

ctrllr.tethers.transformMat = [1 .5 -.5; 1 -.5 0; 1 .5 .5];

ctrllr.tethers.altiTetherKp = 1;
ctrllr.tethers.altiTetherKi = 1;
ctrllr.tethers.altiTetherKd = 1;
ctrllr.tethers.altiTetherTau = 1;

ctrllr.tethers.pitchTetherKp = 1;
ctrllr.tethers.pitchTetherKi = 1;
ctrllr.tethers.pitchTetherKd = 1;
ctrllr.tethers.pitchTetherTau = 1;

ctrllr.tethers.rollTetherKp = 1;
ctrllr.tethers.rollTetherKi = 1;
ctrllr.tethers.rollTetherKd = 1;
ctrllr.tethers.rollTetherTau = 1;

ctrllr.controlSurfaces.aileronKp = 1;
ctrllr.controlSurfaces.aileronKi = 0;
ctrllr.controlSurfaces.aileronKd = 0;
ctrllr.controlSurfaces.aileronTau = 0;
ctrllr.controlSurfaces.aileronMaxDef = 30;


ctrllr.controlSurfaces.elevatorKp = 1;
ctrllr.controlSurfaces.elevatorKi = 0;
ctrllr.controlSurfaces.elevatorKd = 0;
ctrllr.controlSurfaces.elevatorTau = 0;

ctrllr.controlSurfaces.elevatorMaxDef = 30;

%% test signal
alti_sp = 0;
pitch_sp = 0*pi/180;
roll_sp = 1*pi/180;

Rcm_o = zeros(3,1);
Vcm_o = zeros(3,1);
euler = [0;0;0]*pi/180;
OwB = [0;0;0];

flow = [0;0;0];










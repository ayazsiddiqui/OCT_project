clear
clc
format compact
close all

vhcl = PLT.vehicle;

vhcl.setLengthScale(1/2,'');
vhcl.setDensityScale(1,'');
vhcl.setNumTethers(3,'');
vhcl.setNumTurbines(2,'');
vhcl.setBuoyFactor(1.4,'');

% volume and inertias
vhcl.setVolume(0.954,'m^3');
vhcl.setIxx(1,'kg*m^2');
vhcl.setIyy(1,'kg*m^2');
vhcl.setIzz(1,'kg*m^2');
vhcl.setIxy(1,'kg*m^2');
vhcl.setIxz(1,'kg*m^2');
vhcl.setIyz(1,'kg*m^2');
vhcl.setRcb_cm([0;0;0],'m')

% wing
vhcl.setRwingLE_cm([-1.5;0;0],'m');
vhcl.setWingChord(1,'m');
vhcl.setWingAR(10,'');
vhcl.setWingTR(0.8,'');
vhcl.setWingSweep(10,'deg');
vhcl.setWingDihedral(2,'deg');
vhcl.setWingIncidence(0,'deg');
vhcl.setWingNACA('2412','');
vhcl.setWingClMax(1.7,'');
vhcl.setWingClMin(-1.7,'');

% H-stab
vhcl.setRhsLE_wingLE([6;0;0],'m');
vhcl.setHsChord(0.5,'m');
vhcl.setHsAR(8,'');
vhcl.setHsTR(0.8,'');
vhcl.setHsSweep(10,'deg');
vhcl.setHsDihedral(0,'deg');
vhcl.setHsIncidence(0,'deg');
vhcl.setHsNACA('0015','');
vhcl.setHsClMaxl(1.7,'');
vhcl.setHsClMin(-1.7,'');

% V-stab
vhcl.setRvs_wingLE([6;0;0],'m');
vhcl.setVsChord(0.5,'m');
vhcl.setVsSpan(2.5,'m');
vhcl.setVsTR(0.8,'');
vhcl.setVsSweep(10,'deg');
vhcl.setVsNACA('0015','');
vhcl.setVsClMax(1.7,'');
vhcl.setVsClMin(-1.7,'');

% initial conditions
vhcl.setInitialCmPos([0;0;100],'m');
vhcl.setInitialCmVel([0;0;0],'m/s');
vhcl.setInitialEuler([0;0;0],'rad');
vhcl.setInitialAngVel([0;0;0],'rad/s');


vhcl.plot



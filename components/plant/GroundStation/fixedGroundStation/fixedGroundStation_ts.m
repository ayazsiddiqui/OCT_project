clear
clc
format compact

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

%% lifiting body
vhcl = PLT.vehicle;

vhcl.setLengthScale(lengthScale,'');
vhcl.setDensityScale(densityScale,'');
vhcl.setNumTethers(numTethers,'');
vhcl.setNumTurbines(numTurbines,'');
vhcl.setBuoyFactor(1.25,'');

% % % volume and inertias
vhcl.setVolume(945352023.474*1e-9,'m^3');
vhcl.setIxx(6.303080401918E+09*1e-6,'kg*m^2');
vhcl.setIyy(2080666338.077*1e-6,'kg*m^2');
vhcl.setIzz(8.320369733598E+09*1e-6,'kg*m^2');
vhcl.setIxy(0,'kg*m^2');
vhcl.setIxz(81875397.942*1e-6,'kg*m^2');
vhcl.setIyz(0,'kg*m^2');
vhcl.setRcb_cm([0;0;0],'m');

% % % data file name
vhcl.setFluidCoeffsFileName('somefile','');

% % % wing
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

% % % H-stab
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

% % % V-stab
vhcl.setRvs_wingLE([6;0;0],'m');
vhcl.setVsChord(0.5,'m');
vhcl.setVsSpan(2.5,'m');
vhcl.setVsTR(0.8,'');
vhcl.setVsSweep(10,'deg');
vhcl.setVsNACA('0015','');
vhcl.setVsClMax(1.7,'');
vhcl.setVsClMin(-1.7,'');

% % % initial conditions
vhcl.setInitialCmPos([0;0;100],'m');
vhcl.setInitialCmVel([0;0;0],'m/s');
vhcl.setInitialEuler([0;0;0],'rad');
vhcl.setInitialAngVel([0;0;0],'rad/s');

% % % scale the vehicle
vhcl.scaleVehicle

% % % load/generate fluid dynamic data
vhcl.calcFluidDynamicCoefffs


%% ground station
gnd = PLT.gndStn;

gnd.setLengthScale(lengthScale,'');
gnd.setDensityScale(densityScale,'');
gnd.setNumTethers(numTethers,'');

gnd.setIzz(100,'kg*m^2');
gnd.setDampingCoeff(10,'N*m*s');
gnd.setFreeSpinSwitch(0,'');

gnd.setThrAttchPts(vhcl);

% % % initial conditions
gnd.setInitialEuler(0,'rad');
gnd.setInitialAngVel(0,'rad/s');

gnd.scaleGndStn;

%% test signals
f1 = [1;1;0];
f2 = [0;0;0];
f3 = [0;0;0];

attp = gnd.thrAttchPts.Value;

tg1 = cross(attp(:,1),f1);
tg2 = cross(attp(:,2),f2);
tg3 = cross(attp(:,3),f3);

FG = [tg1 tg2 tg3];

test_forces = [f1 f2 f3];

%% run simulink
sim_time = 20;

sim('fixedGroundStation_th')

parseLogsout

tsc.Pos;
tsc.Vel;

n_tet = numTethers;

pos = cell(n_tet,1);

for ii = 1:n_tet
    pos{ii} = squeeze(tsc.Pos.Data(:,ii,:));
    
    for jj = 1:3
        subplot(n_tet,1,ii)
        plot(tsc.Pos.Time,pos{ii}(jj,:));
        if jj == 1
            hold on
            grid on
        elseif jj == 3
            legend('x','y','z');
        end
    end
    
end







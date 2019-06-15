% close
close all;clear;clc;format compact;fclose all;

%% load lookup table
load('dsgnTest_1_lookupTables');

Sref = dsgnData.Sref;
Cref = dsgnData.Cref;
Bref = dsgnData.Bref;


%% parameters
rho_fluid = 1;

%% signals
vel_cm = [0;0;0];
vel_flow = [1;0;0];
euler = [0;0;0];

flpDefl_deg = 0;
ailDefl_deg = 0;

elevDefl_deg = 0;
rudDefl_deg = 0;

%% do a flow sweep
n_step = 21;
alphas = linspace(-25,25,n_step);
betas = linspace(-25,25,n_step);

flow_sweep = zeros(3,n_step);

% switch between alpha or beta test
switch_alp_beta = 2;

for ii = 1:n_step
    if switch_alp_beta == 1
        flow_sweep(:,ii) = [cosd(alphas(ii)) 0 -sind(alphas(ii));0  1 0;...
            sind(alphas(ii))  0 cosd(alphas(ii))]*vel_flow;
    elseif switch_alp_beta == 2
        flow_sweep(:,ii) = [cosd(betas(ii)) sind(betas(ii)) 0;...
            -sind(betas(ii)) cosd(betas(ii)) 0; 0 0 1]*vel_flow;
    end
end
timeVec = 0:length(alphas)-1;

flow_sweep = timeseries(flow_sweep,timeVec);

%% simulate
sim('avlAerodynamics2D_th')
parseLogsout

figure
% Plot things vs alpha
subplot(2,3,1)
plot(squeeze(tsc.angleOfAttackDeg.Data),squeeze(tsc.liftCoeff.Data))
xlabel('$\alpha$, [deg]')
ylabel('Lift Coefficient $C_L$')
grid on

subplot(2,3,2)
plot(squeeze(tsc.angleOfAttackDeg.Data),squeeze(tsc.dragCoeff.Data))
xlabel('$\alpha$, [deg]')
ylabel('Drag Coefficient $C_D$')
grid on


subplot(2,3,3)
plot(squeeze(tsc.angleOfAttackDeg.Data),squeeze(tsc.pitchMomentCoeff.Data))
xlabel('$\alpha$, [deg]')
ylabel({'Pitching Moment','Coefficient $C_m$'})
grid on

subplot(2,3,4)
plot(squeeze(tsc.dragCoeff.Data),squeeze(tsc.liftCoeff.Data))
xlabel('$C_D$')
ylabel('$C_L$')
grid on

subplot(2,3,5)
plot(squeeze(tsc.sideSlipAngleDeg.Data),squeeze(tsc.rollMomentCoeff.Data))
xlabel('$\beta$, [deg]')
ylabel({'Roll Moment','Coefficient $C_l$'})
grid on

subplot(2,3,6)
plot(squeeze(tsc.sideSlipAngleDeg.Data),squeeze(tsc.yawMomentCoeff.Data))
xlabel('$\beta$, [deg]')
ylabel({'Yaw Moment','Coefficient $C_n$'})
grid on


avlPlot_2D_Polars('dsgnTest_1_lookupTables')


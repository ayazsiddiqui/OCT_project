% post processing
% close all

%% colors and linewidth
red = 1/255*[228,26,28];
black = 1/255*[0,0,0];
line_wd = 1;

% % % parse the logged data signals
parseLogsout

%% Scale factors
Lscale = lengthScale;
Dscale = densityScale;

%% resample data
resampleDataRate = 1*Lscale^0.5;
% % % filename = 'testAnimated.gif';
signals = fieldnames(tsc);

newTimeVec = 0:resampleDataRate:tsc.(signals{1}).Time(end);

for ii = 1:length(signals)
    tscResample.(signals{ii}) = resample(tsc.(signals{ii}),newTimeVec);
end

% % % extract the important variables into dummy variables
time = tscResample.inertialCmPos.Time.*(1/Lscale^0.5);
sol_Rcm_o = squeeze(tscResample.inertialCmPos.Data).*(1/Lscale);
sol_Vcmo = squeeze(tscResample.inertialCmVel.Data).*(1/Lscale^0.5);
sol_euler = squeeze(tscResample.eulerAngles.Data);
sol_OwB = squeeze(tscResample.angularVel.Data).*(Lscale^0.5);

%% plot states
plotProps{1} = 'rgb';
if Lscale == 1 && Dscale == 1
    plotProps{2} = '-';
elseif Lscale ~= 1 && Dscale == 1
    plotProps{2} = '--';
elseif Dscale ~= 1 && Lscale == 1
    plotProps{2} = ':';
elseif Dscale ~= 1 && Lscale ~= 1
    plotProps{2} = '.-';
end



ss = get(0,'ScreenSize');
ss = [ss(3) ss(4)];
fig_wid = 1.5*560;
fig_hgt = 1.5*420;
max_horz = floor(ss(1)/fig_wid);
max_vert = floor(ss(2)/fig_hgt);
locs = zeros(max_horz*max_vert,4);
kk = 1;
for jj = 1:max_vert
    for ii = 1:max_horz
        locs(kk,:) = [(ii-1)*fig_wid  ss(2)-(1.2*fig_hgt*jj) fig_wid fig_hgt ];
        kk = kk+1;
    end
end
locs = repmat(locs,5,1);

% % % cm position and set points
fn = 1;
figure(fn);
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_Rcm_o,plotProps,...
    {'$x_{cm}$','$y_{cm}$','$z_{cm}$'},'Position (m)','CM position');
subplot(3,1,3)
% % % setpoint
if numTethers == 1
    plot(time,altiMin*ones(size(time)).*(1/Lscale),'k--',...
    'DisplayName','$Z_{min}$');
    plot(time,altiMax*ones(size(time)).*(1/Lscale),'k--',...
    'DisplayName','$Z_{max}$');
else
plot(time,squeeze(tscResample.altitudeSetpoint.Data.*(1/Lscale)),'k--',...
    'DisplayName','$Z_{sp}$');
end

% % % cm velocity
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_Vcmo,plotProps,...
    {'$V_{x}$','$V_{y}$','$V_{z}$'},'Velocity (m/s)','CM velocity');

% % % euler angles
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_euler*180/pi,plotProps,...
    {'$\phi$','$\theta$','$\psi$'},'Angle (deg)','Euler Angles');
% % % setpoints
subplot(3,1,1)
plot(time,squeeze(tscResample.rollSetpoint.Data)*180/pi,'k--',...
    'DisplayName','$\phi_{sp}$');
subplot(3,1,2)
plot(time,squeeze(tscResample.pitchSetpoint.Data)*180/pi,'k--',...
    'DisplayName','$\theta_{sp}$');


% % % angular velocities
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,sol_OwB,plotProps,...
%     {'$\omega_{x}$','$\omega_{y}$','$\omega_{z}$'},'Ang vel (rad/s)','Angular velocities');

%% other angles
% elevation angle
elevAngle = (180/pi)*atan2(sol_Rcm_o(3,:),sqrt(sum(sol_Rcm_o(1:2,:).^2,1)));
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))

% azimuth angle
azimuthAngle = (180/pi)*atan2(sol_Rcm_o(2,:),sol_Rcm_o(1,:));

vectorPlotter(time,[elevAngle;azimuthAngle],plotProps,...
    {'Elevation','Azimuth'},'Angle (deg)','Other angles');







%% plot forces
% % % % buoyancy force
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyBuoyForce.Data.*(1/Lscale^3),plotProps,...
%     {'$F_{buoy,x}$','$F_{buoy,y}$','$F_{buoy,z}$'},'Force (N)','Buoyancy Force');
% 
% % % % gravity force
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyGravForce.Data.*(1/Lscale^3),plotProps,...
%     {'$F_{grav,x}$','$F_{grav,y}$','$F_{grav,z}$'},'Force (N)','Gravitational Force');
% 
% % % % Aero force
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyAeroForce.Data.*(1/Lscale^3),plotProps,...
%     {'$F_{aero,x}$','$F_{aero,y}$','$F_{aero,z}$'},'Force (N)','Aero Force');
% 
% % % % total force
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyTotForce.Data.*(1/Lscale^3),plotProps,...
%     {'$F_{tot,x}$','$F_{tot,y}$','$F_{tot,z}$'},'Force (N)','Total Force');



%% plot moments
% % % % aero moment
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyAeroMoment.Data.*(1/Lscale^4),plotProps,...
%     {'$M_{aero,x}$','$M_{aero,y}$','$M_{aero,z}$'},'Moment (N-m)','Aero Moment');
% 
% % % % tether moment
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyTetherMoment.Data.*(1/Lscale^4),plotProps,...
%     {'$M_{tether,x}$','$M_{tether,y}$','$M_{tether,z}$'},'Moment (N-m)','Tether Moment');
% 
% % % % total moment
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyTotMoment.Data.*(1/Lscale^4),plotProps,...
%     {'$M_{total,x}$','$M_{total,y}$','$M_{total,z}$'},'Moment (N-m)','Total Moment');

% % % % turbine moment
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyTurbineMoment.Data.*(1/Lscale^4),plotProps,...
%     {'$M_{turb,x}$','$M_{turb,y}$','$M_{turb,z}$'},'Moment (N-m)','Turbine Moment');
% 
% % % % buoyancy moment
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyBuoyMoment.Data.*(1/Lscale^4),plotProps,...
%     {'$M_{buoy,x}$','$M_{buoy,y}$','$M_{buoy,z}$'},'Moment (N-m)','Buoyancy Moment');

%% plot control signals
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.thrReleseSpeeds.Data.*(1/Lscale^0.5),plotProps,...
    {'$u_{port}$','$u_{aft}$','$u_{stbd}$'},'Speed (m/s)','Tether release speeds');

fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.thrLengths.Data.*(1/Lscale),plotProps,...
    {'$L_{port}$','$L_{aft}$','$L_{stbd}$'},'Length (m)','Tether lengths');


fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.Csdeflectns.Data,plotProps,...
    {'$\delta_{port-alrn}$','$\delta_{stbd-alrn}$','$\delta_{elevator}$','$\delta_{rudder}$'},...
    'Angle (deg)','Control surface defelctions');



%% local forces
% angle of attack
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,squeeze(tscResample.angleOfAtk.Data),plotProps,...
    {'Port wing','Stbd wing','H-stab','V-stab'},...
    'Angle (deg)','Angle of attack');

% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,squeeze(tscResample.liftCoeff.Data),plotProps,...
%     {'Port wing','Stbd wing','H-stab','V-stab'},...
%     'CL','Lift coefficient');
% 
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,squeeze(tscResample.dragCoeff.Data),plotProps,...
%     {'Port wing','Stbd wing','H-stab','V-stab'},...
%     'CD','Drag coefficient');

% 
% surfNames = {'Port wing','Stbd wing','H-stab','V-stab'};
% for ii = 1:4
%     fn = fn+1;
%     figure(fn)
%     set(gcf,'Position',locs(fn,:))
%     vectorPlotter(time,tscResample.bdyLiftForce.Data(:,ii,:),plotProps,...
%         {'$F_{x}$','$F_{y}$','$F_{z}$'},'F(N)',strcat(surfNames{ii},' Lift'));
% 
% end
    



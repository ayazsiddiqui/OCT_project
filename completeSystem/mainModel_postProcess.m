% post processing
% close all

%% colors and linewidth
red = 1/255*[228,26,28];
black = 1/255*[0,0,0];
line_wd = 1;

% parse the logged data signals
parseLogsout

%% resample data
resampleDataRate = 0.1;
% filename = 'testAnimated.gif';
signals = fieldnames(tsc);

newTimeVec = 0:resampleDataRate:tsc.(signals{1}).Time(end);

for ii = 1:length(signals)
    tscResample.(signals{ii}) = resample(tsc.(signals{ii}),newTimeVec);
end

time = tscResample.inertialCmPos.Time;

sol_Ri_o = tscResample.allNodePos.Data;
sol_Rcm_o = tscResample.inertialCmPos.Data;
sol_Vcmo = tscResample.inertialCmVel.Data;
sol_euler = tscResample.eulerAngles.Data;
sol_OwB = tscResample.angularVel.Data;

%% plot states
plotProps{1} = 'rgb';
if run_no == 1
    plotProps{2} = '-';
elseif run_no == 2
    plotProps{2} = '--';
end

ss = get(0,'ScreenSize');
ss = [ss(3) ss(4)];
fig_wid = 560;
fig_hgt = 420;
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

% cm position and set points
fn = 1;
figure(fn);
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_Rcm_o,plotProps,...
    {'$x_{cm}$','$y_{cm}$','$z_{cm}$'},'CM position (m)');
subplot(3,1,3)
plot(time,squeeze(tscResample.altitudeSetpoint.Data),'k--',...
    'DisplayName','$Z_{sp}$');

% cm velocity
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_Vcmo,plotProps,...
    {'$V_{x}$','$V_{y}$','$V_{z}$'},'CM velocity (m/s)');

% euler angles
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_euler*180/pi,plotProps,...
    {'$\phi$','$\theta$','$\psi$'},'Euler Ang (deg)');
subplot(3,1,1)
plot(time,squeeze(tscResample.rollSetpoint.Data)*180/pi,'k--',...
    'DisplayName','$\phi_{sp}$');
subplot(3,1,2)
plot(time,squeeze(tscResample.pitchSetpoint.Data)*180/pi,'k--',...
    'DisplayName','$\theta_{sp}$');


% angular velocities
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,sol_OwB,plotProps,...
    {'$\omega_{x}$','$\omega_{y}$','$\omega_{z}$'},'Ang vel (rad/s)');

%% plot forces
% buoyancy force
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyBuoyForce.Data,plotProps,...
    {'$F_{buoy,x}$','$F_{buoy,y}$','$F_{buoy,z}$'},'Force (N)');

% gravity force
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyGravForce.Data,plotProps,...
    {'$F_{grav,x}$','$F_{grav,y}$','$F_{grav,z}$'},'Force (N)');

% gravity force
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyAeroForce.Data,plotProps,...
    {'$F_{aero,x}$','$F_{aero,y}$','$F_{aero,z}$'},'Force (N)');

% gravity force
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyTotForce.Data,plotProps,...
    {'$F_{tot,x}$','$F_{tot,y}$','$F_{tot,z}$'},'Force (N)');



%% plot moments
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyAeroMoment.Data,plotProps,...
    {'$M_{aero,x}$','$M_{aero,y}$','$M_{aero,z}$'},'Aero Moment (N-m)');

fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyTetherMoment.Data,plotProps,...
    {'$M_{tether,x}$','$M_{tether,y}$','$M_{tether,z}$'},'Tether Moment (N-m)');

fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.bdyTotMoment.Data,plotProps,...
    {'$M_{total,x}$','$M_{total,y}$','$M_{total,z}$'},'Total Moment (N-m)');



% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyTurbineMoment.Data,plotProps,...
%     {'$M_{turb,x}$','$M_{turb,y}$','$M_{turb,z}$'},'Turbine Moment (N-m)');
% 
% fn = fn+1;
% figure(fn)
% set(gcf,'Position',locs(fn,:))
% vectorPlotter(time,tscResample.bdyBuoyMoment.Data,plotProps,...
%     {'$M_{buoy,x}$','$M_{buoy,y}$','$M_{buoy,z}$'},'Buoyancy Moment (N-m)');

%% plot control signals
fn = fn+1;
figure(fn)
set(gcf,'Position',locs(fn,:))
vectorPlotter(time,tscResample.thrReleseSpeeds.Data,plotProps,...
    {'$u_{port}$','$u_{aft}$','$u_{stbd}$'},'Tether release speeds (m/s)');



%% animations plots

nNodes = tp.tethers(1).numNodes;
nTethers = length(tp.tethers);

s_R = cell(nTethers,1);
s_Rn_o = cell(nTethers,1);
s_R1_o = cell(nTethers,1);

for ii = 1:nTethers
    s_R{ii} = squeeze(sol_Ri_o(:,:,ii,:));
    s_R1_o{ii} = s_R{ii}(:,1,:);
    s_Rn_o{ii} = s_R{ii}(:,end,:);
end

bx = zeros(nTethers,2);
by = zeros(nTethers,2);
bz = zeros(nTethers,2);

for ii = 1:nTethers
    [xmin,xmax] = bounds(squeeze(s_R{ii}(1,:,:)),'all');
    [ymin,ymax] = bounds(squeeze(s_R{ii}(2,:,:)),'all');
    [zmin,zmax] = bounds(squeeze(s_R{ii}(3,:,:)),'all');
    
    bx(ii,:) = round([floor(xmin-5),ceil(xmax+5)],-1);
    by(ii,:) = round([floor(ymin-5),ceil(ymax+5)],-1);
    bz(ii,:) = round([floor(zmin-5),ceil(zmax+5)],-1);
end

%% plot
n_steps = length(time);
if plot_animation == 0
    return
end
fn = fn+1;
figure(fn)

% video setting
video = VideoWriter('vid_Test', 'Motion JPEG AVI');
video.FrameRate = 1/resampleDataRate;

mov(1:n_steps)=struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren');

for ii = 1:n_steps
    
    if ii > 1
        h = findall(gca,'type','line','color',red,'-or','color',black);
        delete(h);
    end
    
    for kk = 1:nTethers
        p3d_1 = plot3(s_R{kk}(1,:,ii),s_R{kk}(2,:,ii),s_R{kk}(3,:,ii),...
            '-+','linewidth',line_wd,'color',black);
        hold on
        pRcm_n = plot3([s_R{kk}(1,end,ii) sol_Rcm_o(1,1,ii)],...
            [s_R{kk}(2,end,ii) sol_Rcm_o(2,1,ii)],...
            [s_R{kk}(3,end,ii) sol_Rcm_o(3,1,ii)],...
            '-','linewidth',line_wd,'color',red);
    end
    
    if ii == 1
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
        xlim([-max(abs(bx(:))) max(abs(bx(:)))]);
        ylim([-max(abs(by(:))) max(abs(by(:)))]);
        zlim([0 max(bz(:))]);
%         axis equal
        hold on
        grid on
    end
    
    title(['Time = ',sprintf('%0.2f', time(ii)),' s'])
    
    try
%         waitforbuttonpress
    catch
        break
    end
    F(ii) = getframe(gcf);

end

if make_video == 1
    open(video)
    for i = 1:length(F)
        writeVideo(video, F(i));
    end
    close(video)
end
clear
clc
format compact
% close all

% % merged

%% common parameters
lengthScale = 1;
densityScale = 1;
numTethers = 3;
numTurbines = 2;

%% tethers
thr = PLT.tether;s

thr.setLengthScale(lengthScale,'');
thr.setDensityScale(densityScale,'');
thr.setNumTethers(numTethers,'');

thr.setNumNodes(10,'');
thr.setThrDiameter([0.01 0.02 0.01],'m');
thr.setThrDensity(1300*ones(1,numTethers),'kg/m^3');
thr.setThrYoungs(3.8e9*ones(1,numTethers),'N/m^2');
thr.setThrDampingRatio(0.05*ones(1,numTethers),'');
thr.setThrDragCoeff(0.5*ones(1,numTethers),'');

thr.scaleTether;

%% dummy constants
mass = 1000;
rho_fluid = 1000;
grav = 9.81;


%% test signals
dt = 1/100;

amp = 10;
omega = 5;
sim_time = 10;

z_pos = 100;

ini_length = sqrt(z_pos^2 + 2*amp^2);
flow = [0;0;0];

tVec = 0:dt/10:sim_time;

Rn_o_test = timeseries();
Vn_o_test = timeseries();

Rn_o_test.Time = tVec;
Vn_o_test.Time = tVec;

for ii = 1:numTethers
    Rn_o_test.Data(:,ii,:) = [amp*cos(omega*tVec)+ 5*(ii-1);
        amp*sin(omega*tVec);
        z_pos*ones(size(tVec))];
    
    Vn_o_test.Data(:,ii,:) = [-amp*omega*sin(omega*tVec);
        amp*omega*cos(omega*tVec);
        0*ones(size(tVec))];
    
end


%% set time step and simulate
% simulate
open_system('tetherNnode_th')
sim('tetherNnode_th')




%% post process

% parse data
parseLogsout

% colors
red = 1/255*[228,26,28];
blue = 1/255*[55,126,184];
green = 1/255*[77,175,74];
purple = 1/255*[152,78,163];
line_wd = 0.75;

time = tsc.F1_tet.Time;

sol_Ri_o = tsc.Ri_o.Data;
sol_F1_tet = tsc.F1_tet.Data;
sol_Fn_tet = tsc.Fn_tet.Data;
sol_Vi_o = tsc.Vi_o.Data;

nNodes = thr.numNodes.Value;

s_R = cell(numTethers,1);
s_Rn_o = cell(numTethers,1);
s_R1_o = cell(numTethers,1);

for ii = 1:numTethers
    s_R{ii} = squeeze(sol_Ri_o(:,:,ii,:));
    s_R1_o{ii} = s_R{ii}(:,1,:);
    s_Rn_o{ii} = s_R{ii}(:,end,:);
end

bx = zeros(numTethers,2);
by = zeros(numTethers,2);
bz = zeros(numTethers,2);

for ii = 1:numTethers
[xmin,xmax] = bounds(squeeze(s_R{ii}(1,:,:)),'all'); 
[ymin,ymax] = bounds(squeeze(s_R{ii}(2,:,:)),'all'); 
[zmin,zmax] = bounds(squeeze(s_R{ii}(3,:,:)),'all');

bx(ii,:) = [xmin,xmax];
by(ii,:) = [ymin,ymax];
bz(ii,:) = [zmin,zmax];

end
    

%% make movie
movie_frame_rate = 50;
skip_step = (1/dt)/movie_frame_rate;
t_steps = length(time);

t_f = time(end);

n_step_idx = 1:skip_step:t_steps;
n_steps = length(n_step_idx);
t_snap = time(n_step_idx);

% video setting
video = VideoWriter('vid_Test', 'Motion JPEG AVI');
video.FrameRate = movie_frame_rate;
num_frames = n_steps;

mov(1:n_steps)=struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren')


for ii = 1:n_steps
    
    figure(1)
    if ii > 1
        h = findall(gca, 'type', 'line','color',red);
        delete(h);
    end
    
    for kk = 1:numTethers
        
        p3d_1 = plot3(s_R{kk}(1,:,ii),s_R{kk}(2,:,ii),s_R{kk}(3,:,ii),'-+','color',red);
        hold on
        
    end
    
   if ii == 1
        xlabel('Y (m)'); ylabel('Y (m)'); zlabel('Z (m)')
        
        xlim([min(bx(:)) max(bx(:))]);
        ylim([min(by(:)) max(by(:))]);
        zlim([min(bz(:)) max(bz(:))]);
        hold on
        grid on
   end
    
    title(['Time = ',sprintf('%0.2f', t_snap(ii)),' s'])
    F(ii) = getframe(gcf);
end


% open(video)
% for i = 1:length(F)
%     writeVideo(video, F(i));
% end
% close(video)



%%
% vssBlk ='kelvinVoigtTether_cl/kelvinVoigtTether';





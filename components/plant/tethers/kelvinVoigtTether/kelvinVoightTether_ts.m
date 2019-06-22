clear
clc
format compact

% start
%% environment param
G = 9.81;
rho_f = 1000;

mass = 1000;

%% tether parameters
class_thr(1).numNodes = 2;
class_thr(1).diameter = 0.05;
class_thr(1).youngsModulus = 3.8e9;
class_thr(1).dampingRatio = 0.05;
class_thr(1).dragCoeff = 0.5;
class_thr(1).density = 1300;
class_thr(1).vehicleMass = 100;

class_thr(2).numNodes = 2;
class_thr(2).diameter = 0.05;
class_thr(2).youngsModulus = 3.8e9;
class_thr(2).dampingRatio = 0.05;
class_thr(2).dragCoeff = 0.5;
class_thr(2).density = 1300;
class_thr(2).vehicleMass = 100;

% class_thr(3).numNodes = 10;
% class_thr(3).diameter = 0.05;
% class_thr(3).youngsModulus = 3.8e9;
% class_thr(3).dampingRatio = 0.05;
% class_thr(3).dragCoeff = 0.5;
% class_thr(3).density = 1300;
% class_thr(3).vehicleMass = 100;


%% test signals
dt = 1/1000;

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

for ii = 1:length(class_thr)
    Rn_o_test.Data(:,ii,:) = [amp*cos(omega*tVec)+ 5*(ii-1);
        amp*sin(omega*tVec);
        z_pos*ones(size(tVec))];
    
    Vn_o_test.Data(:,ii,:) = [-amp*omega*sin(omega*tVec);
        amp*omega*cos(omega*tVec);
        0*ones(size(tVec))];
    
end


%% set time step and simulate
% simulate

open_system('kelvinVoigtTether_th')

kelvinVoightTether_init

try
    sim('kelvinVoigtTether_th')
catch
    pause(3)
    sim('kelvinVoigtTether_th')
end

save_system('kelvinVoigtTether_th')
close_system('kelvinVoigtTether_th')



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

nNodes = class_thr(1).numNodes;
nTethers = length(class_thr);

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
    
    for kk = 1:nTethers
        
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





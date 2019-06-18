clear
clc
format compact

% start
%% environment param
G = 9.81;
rho_f = 1000;

mass = 1000;

amp = 10;
omega = 1;
sim_time = 5;


%% tether parameters
class_thr(1).R1_g = [0; -1; 0];
class_thr(1).Rn_cm = [0; -1; -1];
class_thr(1).numNodes = 10;
class_thr(1).diameter = 0.05;
class_thr(1).youngsModulus = 3.8e9;
class_thr(1).dampingRatio = 0.05;
class_thr(1).dragCoeff = 0.5;
class_thr(1).density = 1300;
class_thr(1).vehicleMass = 100;
class_thr(1).ini_R1_o = [0; 0; 0];
class_thr(1).ini_Rn_o = [0; 0; 100];


class_thr(2).R1_g = [1; 0; 0];
class_thr(2).Rn_cm = [1; 0; -1];
class_thr(2).numNodes = 10;
class_thr(2).diameter = 0.05;
class_thr(2).youngsModulus = 3.8e9;
class_thr(2).dampingRatio = 0.05;
class_thr(2).dragCoeff = 0.5;
class_thr(2).density = 1300;
class_thr(2).vehicleMass = 100;
class_thr(2).ini_R1_o = [1; 0; 0];
class_thr(2).ini_Rn_o = [1; 0; 100];
%
% class_thr(3).R1_g = [0; 1; 0];
% class_thr(3).Rn_cm = [0; 1; -1];
% class_thr(3).numNodes = 2;
% class_thr(3).diameter = 0.05;
% class_thr(3).youngsModulus = 3.8e9;
% class_thr(3).dampingRatio = 0.05;
% class_thr(3).dragCoeff = 0.5;
% class_thr(3).density = 1300;
% class_thr(3).vehicleMass = 100;
% class_thr(3).ini_R1_o = [0; 1; 0];
% class_thr(3).ini_Rn_o = [0; 1; 100];

for ii = 1:length(class_thr)
    
    class_thr(ii).initNodePoss = [...
        linspace(class_thr(ii).ini_R1_o(1),class_thr(ii).ini_Rn_o(1),class_thr(ii).numNodes);...
        linspace(class_thr(ii).ini_R1_o(2),class_thr(ii).ini_Rn_o(2),class_thr(ii).numNodes);...
        linspace(class_thr(ii).ini_R1_o(3),class_thr(ii).ini_Rn_o(3),class_thr(ii).numNodes)];
    
    class_thr(ii).initNodePoss = class_thr(ii).initNodePoss(:,2:end-1);
    class_thr(ii).initNodeVels = zeros(size(class_thr(ii).initNodePoss));
    
    n_ini_pos(:,:,ii) = class_thr(ii).initNodePoss;
    n_ini_vel(:,:,ii) = class_thr(ii).initNodeVels;
    
end

%% set time step and simulate
dt = 1/1000;
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

s_R = cell(nNodes,nTethers);

for kk = 0:nTethers-1
    for jj = 1:nNodes
        for ii = 1:length(time)
            s_R_int(ii,:) =  sol_Ri_o((3*kk+1):(3*kk+3),jj,ii)';
        end
        s_R{jj,kk+1} =  s_R_int;
    end
end

s_Rn_o = s_R(end,:);
s_R1_o = s_R(1,:);

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

% separate x,y and z cordinates
p3x = NaN(nNodes,n_steps,nTethers);
p3y = NaN(nNodes,n_steps,nTethers);
p3z = NaN(nNodes,n_steps,nTethers);

for kk = 1:nTethers
    for jj = 1:nNodes
        int_p1 = s_R{jj,kk};
        p3x(jj,:,kk) = int_p1(n_step_idx,1)';
        p3y(jj,:,kk) = int_p1(n_step_idx,2)';
        p3z(jj,:,kk) = int_p1(n_step_idx,3)';
    end
end

[Sx,Lx] = bounds(p3x,'all');
Sx = Sx(1); Lx = Lx(1);
[Sy,Ly] = bounds(p3y,'all');
Sy = Sy(1)-1; Ly = Ly(1)+1;
[Sz,Lz] = bounds(p3z,'all');
Sz = Sz(1); Lz = Lz(1);

% axis square
% axis equal

for ii = 1:n_steps
    
    figure(1)
    if ii > 1
        h = findall(gca, 'type', 'line','color',red);
        delete(h);
    end
    
    for kk = 1:nTethers
        
        p3d_1 = plot3(p3x(:,ii,kk),p3y(:,ii,kk),p3z(:,ii,kk),'-+','color',red);
        hold on
        
    end
    
   if ii == 1
        xlabel('Y (m)'); ylabel('Y (m)'); zlabel('Z (m)')
        
        xlim([Sx Lx]); ylim([Sy Ly]); zlim([Sz Lz]);
        hold on
        grid on
   end
    
    title(['Time = ',sprintf('%0.2f', t_snap(ii)),' s'])
    
    F(ii) = getframe(gcf);
end

%%
% vssBlk ='kelvinVoigtTether_cl/kelvinVoigtTether';





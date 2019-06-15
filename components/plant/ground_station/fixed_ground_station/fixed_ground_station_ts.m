clear
clc
format compact

%% mask parameters
Izz = 1;
c_damp = 1;

class_thr(1).R1_g = [0; -1; 0];
class_thr(1).Rn_cm = [0; -0.5; -1];

% class_thr(2).R1_g = [1; 0; 0];
% class_thr(2).Rn_cm = [1; 0; -1];
% 
% class_thr(3).R1_g = [0; 1; 0];
% class_thr(3).Rn_cm = [0; 0.5; -1];

% initial conditions
ini_OwP = 0;
ini_psiPlat = 0*pi/180;

f1 = [1;1;0];
f2 = [0;0;0];
f3 = [0;0;0];


tg1 = cross(class_thr(1).R1_g,f1);
tg2 = cross(class_thr(2).R1_g,f2);
tg3 = cross(class_thr(3).R1_g,f3);

tf1 = cross(class_thr(1).Rn_cm,f1);
tf2 = cross(class_thr(2).Rn_cm,f2);
tf3 = cross(class_thr(3).Rn_cm,f3);

FG = [tg1 tg2 tg3];


TF = [tf1 tf2 tf3];


test_forces = [f1 f2 f3];

%% run simulink
sim_time = 20;

sim('fixed_ground_station_th')

parseLogsout

tsc.F_out;
tsc.M_out;
tsc.Pos;
tsc.Vel;

n_tet = length(class_thr);

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







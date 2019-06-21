clear
clc
format compact

%% mask parameters
Izz = 1;
c_damp = 1;

class_gndStn(1).tetherAttchPt = [0; -1; 0];
class_gndStn(2).tetherAttchPt = [1; 0; 0];
class_gndStn(3).tetherAttchPt = [0; 1; 0];

% initial conditions
ini_OwP = 0;
ini_psiPlat = 0*pi/180;

f1 = [1;1;0];
f2 = [0;0;0];
f3 = [0;0;0];


tg1 = cross(class_gndStn(1).tetherAttchPt,f1);
tg2 = cross(class_gndStn(2).tetherAttchPt,f2);
tg3 = cross(class_gndStn(3).tetherAttchPt,f3);

FG = [tg1 tg2 tg3];

test_forces = [f1 f2 f3];

%% run simulink
sim_time = 20;

sim('fixedGroundStation_th')

parseLogsout

tsc.Pos;
tsc.Vel;

n_tet = length(class_gndStn);

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







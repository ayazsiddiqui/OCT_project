clear
clc
format compact

% start
%% environment param
G = 9.81;
rho_f = 1000;

mass = 1000;




%% tether parameters
class_thr(1).R1_g = [0; -1; 0];
class_thr(1).Rn_cm = [0; -1; -1];
class_thr(1).numNodes = 4;
class_thr(1).diameter = 0.05;
class_thr(1).youngsModulus = 3.8e9;
class_thr(1).dampingRatio = 0.05;
class_thr(1).dragCoeff = 0.5;
class_thr(1).density = 1300;
class_thr(1).vehicleMass = 100;
class_thr(1).ini_R1_o = [0; 0; 0];
class_thr(1).ini_Rn_o = [0; 0; 100];

% 
% class_thr(2).R1_g = [1; 0; 0];
% class_thr(2).Rn_cm = [1; 0; -1];
% class_thr(2).numNodes = 2;
% class_thr(2).diameter = 0.05;
% class_thr(2).youngsModulus = 3.8e9;
% class_thr(2).dampingRatio = 0.05;
% class_thr(2).dragCoeff = 0.5;
% class_thr(2).density = 1300;
% class_thr(2).vehicleMass = 100;
% class_thr(2).ini_R1_o = [1; 0; 0];
% class_thr(2).ini_Rn_o = [1; 0; 100];
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

%% simulate
sim('kelvinVoigtTether_th2')

%% post process
parseLogsout

time = tsc.F1_tet.Time;
s_F1_tet = tsc.F1_tet.Data;
s_Fn_tet = tsc.Fn_tet.Data;
s_Ri_o = tsc.Ri_o.Data;
s_Vi_o = tsc.Vi_o.Data;



%%
% vssBlk ='kelvinVoigtTether_cl/kelvinVoigtTether';
% oldVariant = get_param(vssBlk,'OverrideUsingVariant');
% 
% if any([class_thr.numNodes] == 2)
%     set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNode2_variant');
% else
%     set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNodeN_variant');
% end
% 
% sim('kelvinVoigtTether_th')




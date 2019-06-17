clear
clc
format compact

%% start
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
class_thr(1).ini_R1_o = [0; -1; 0];
class_thr(1).ini_Rn_o = [0; -1; 100];


class_thr(2).R1_g = [1; 0; 0];
class_thr(2).Rn_cm = [1; 0; -1];
class_thr(2).numNodes = 2;
class_thr(2).diameter = 0.05;
class_thr(2).youngsModulus = 3.8e9;
class_thr(2).dampingRatio = 0.05;
class_thr(2).dragCoeff = 0.5;
class_thr(2).density = 1300;
class_thr(2).vehicleMass = 100;
class_thr(2).ini_R1_o = [1; 0; 0];
class_thr(2).ini_Rn_o = [1; 0; 100];

class_thr(3).R1_g = [0; 1; 0];
class_thr(3).Rn_cm = [0; 1; -1];
class_thr(3).numNodes = 4;
class_thr(3).diameter = 0.05;
class_thr(3).youngsModulus = 3.8e9;
class_thr(3).dampingRatio = 0.05;
class_thr(3).dragCoeff = 0.5;
class_thr(3).density = 1300;
class_thr(3).vehicleMass = 100;
class_thr(3).ini_R1_o = [0; 1; 0];
class_thr(3).ini_Rn_o = [0; 1; 100];


%%
numTethers = length(class_thr);

% initialize variables
numNodes = zeros(1,numTethers);

for ii = 1:numTethers
    numNodes(ii) = class_thr(ii).numNodes;
end

numNodes = min(numNodes);


vssBlk ='kelvinVoigtTether_th/kelvinVoigtTether';
oldVariant = get_param(vssBlk,'OverrideUsingVariant');

if any(numNodes==2)
    set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNode2_variant');
else
    set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNodeN_variant');
end

sim('kelvinVoigtTether_th')




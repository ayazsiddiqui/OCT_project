 function [F_tension,magTension]  = fcn2(R1_o,V1_o,Rn_o,Vn_o,Ri_o,Vi_o,thrLength,F_tension,magTension,...
     numNode, massBdy, zetaTether, youngsTether, diaTether)

if numNode == 2
    Ri_o = [R1_o Rn_o];
    Vi_o = [V1_o Vn_o];
else
    Ri_o = [R1_o Ri_o Rn_o];
    Vi_o = [V1_o Vi_o Vn_o];
end
 
% cross sectional area
CS = (pi/4)*diaTether^2;

% total spring stiffness
k_tot = youngsTether*CS/thrLength;

% total damping coefficient
c_tot = zetaTether*2*sqrt(k_tot*massBdy);

% inidividual spring stiffness and damping
kElem = k_tot*(numNode-1);
cElem = c_tot*(numNode-1);
LiElem = thrLength/(numNode-1);

% distance between i node to i+1 node
Relem_o = Ri_o(:,2:end) - Ri_o(:,1:end-1);
distNodes = sqrt(sum(Relem_o.^2,1));

% direction of the tensile force
tensionDirection = Relem_o./repmat(distNodes,3,1);

% velocity of the i^th wrt (i+1)^th node
Velem_o = Vi_o(:,2:end) - Vi_o(:,1:end-1);

% rate of change of distance between nodes
d_dt_distNodes = dot(Relem_o,Velem_o)./distNodes;

% find the links that are in tension
inTension = distNodes > LiElem;

% calculate the tension force
magSpring = kElem*(inTension.*(distNodes - LiElem));
magDamper = cElem*(inTension.*d_dt_distNodes);

Ftensile = (magSpring + magDamper).*tensionDirection;

% distribute the tensile forces among nodes
nodalTension1 = Ftensile(:,1);
nodalTensionInt = Ftensile(:,1:end-1) - Ftensile(:,2:end);
nodalTensionN = -Ftensile(:,end);

% tension forces on each node
F_tension(:,1) = nodalTension1;
F_tension(:,2:end-1) = nodalTensionInt;
F_tension(:,end) = nodalTensionN;

magTension(1) = sqrt(sum(Ftensile(:,1).^2,1));
magTension(2:end-1) = sqrt(sum(Ftensile(:,2:end-1).^2,1));
magTension(end) = sqrt(sum(Ftensile(:,end).^2,1));


end








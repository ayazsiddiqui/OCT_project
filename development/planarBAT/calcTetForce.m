function [sumF,nodeMass] = calcTetForce(Vflow,nodePos,nodeVel,tetLength,numNode, tetDia, tetDensity,tetCD, tetYoungs, tetZeta, mT, rhoF, gravAcc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% downward force
nodeMass = (tetDensity*(pi/4)*tetLength*tetDia^2)/(numNode - 2);
nodeBuoy = (rhoF*(pi/4)*tetLength*tetDia^2)/(numNode - 2);

fGrav = (nodeMass-nodeBuoy)*repmat([0;-gravAcc],1,numNode);

% R_i+1 - Ri
Rip1_i = nodePos(:,2:end) - nodePos(:,1:end-1);
Vip1_i = nodeVel(:,2:end) - nodeVel(:,1:end-1);

% distance between nodes and its derivative
di =  vecnorm(Rip1_i);
diDot = dot(Rip1_i,Vip1_i)./di;

% direction of tensilde force and link angles
direction = Rip1_i./repmat(di,2,1);
angles = atan2(Rip1_i(2,:),Rip1_i(1,:));

% drag force
Vrel = repmat(Vflow,1,numNode-1) - (nodeVel(:,2:end) + nodeVel(:,1:end-1))/2;
mVrel = vecnorm(Vrel);
mVrel(mVrel<eps) = eps;
dragDir = Vrel./repmat(mVrel,2,1);
fDrag = tetCD*0.5*rhoF*tetDia*repmat((di.*(mVrel.^2).*sin(angles)),2,1).*dragDir;
fDrag(isnan(fDrag)) = 0;

% tether forces
kTot = (tetYoungs*(pi/4)*tetDia^2)/tetLength;
ctot = 2*tetZeta*sqrt(kTot*mT);

kLink = kTot*(numNode-1);
cLink = ctot*(numNode-1);

%
Li = tetLength/(numNode-1);
tetF = zeros(size(di));
tester = di>Li;
springF = kLink*(di-Li);
dampF = cLink*diDot;
tetF(tester) = springF(tester) + dampF(tester);

% nodal forces
nodeDrag = (cat(2,[0;0],fDrag) + cat(2,fDrag,[0;0]))/2;
nodeTen = cat(2,repmat(tetF,2,1).*direction,[0;0]) - cat(2,[0;0],repmat(tetF,2,1).*direction);

sumF = fGrav + nodeDrag + nodeTen;


end


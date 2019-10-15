function [nodeForceVecs,gravMass] = nodeForce(flowVel,nodePositions,nodeVelocities,...
    unstretchedLength,mass,zeta,diameter,youngsMod,Cd,rhoF,rhoTet,gravAcc)

% Number of nodes
N = size(nodePositions,2);

% weigth force
totalVolume = (pi/4)*(diameter^2)*unstretchedLength;
gravMass = rhoTet*totalVolume/(N-2);
buoyMass = rhoF*totalVolume/(N-2);
vertForce = (buoyMass-gravMass)*repmat([0;gravAcc],[1 N]);

% total spring stiffness
totalSpringConst = youngsMod*(pi/4)*(diameter^2)/unstretchedLength;

% total damping coefficient
totalDampCoeff = zeta*2*sqrt(totalSpringConst*mass);

% Spring constant and damping coefficient for each link
linkSpringConst  = totalSpringConst*(N-1);
linkDampingCoeff = totalDampCoeff*(N-1);
meanLinkLength   = unstretchedLength/(N-1);

% Vector from one node to another
linkVecs  = diff(nodePositions,1,2);

% angle between linkVecs
thetas = atan2(linkVecs(1,:),linkVecs(2,:));

% Length of each link
linkLength = sqrt(sum(linkVecs.^2));

% projected area
projArea = diameter.*linkLength.*sin(thetas);
linkVel = (nodeVelocities(:,2:end)-nodeVelocities(:,1:end-1))/2;

% apparent flow on link
linkAppFlow = repmat(flowVel,[1 N-1]) - linkVel;

% Link center apparent flow magnitude 3xN-1
linkAppFlowMag = sqrt(sum(linkAppFlow.^2)); % Sum over columns

% dynamic pressure
dynPress = 0.5*rhoF*linkAppFlowMag.^2;

% drag force
linkDragForce = repmat(Cd.*dynPress.*projArea,[2 1]).*(linkAppFlow./repmat(linkAppFlowMag,[2 1]));
% Drag force on nodes
nodeDragForce = [linkDragForce(:,1) (linkDragForce(:,1:end-1)+linkDragForce(:,2:end)) linkDragForce(:,end)]/2;

% Link direction unit vectors
linkUnitVecs = linkVecs./repmat(sqrt(sum(linkVecs.^2,1)),[2 1]);

% Rate of change of vector from one node to another
linkLengthDeriv = dot(nodeVelocities(:,2:end)-nodeVelocities(:,1:end-1),linkUnitVecs);

% Magnitude of spring force on each link
springForces = zeros(size(linkLength));
springForces(linkLength>meanLinkLength) = linkSpringConst*(linkLength(linkLength>meanLinkLength)-meanLinkLength);

% Magnitude of damping force on each link
damperForces  = zeros(size(linkLengthDeriv));
damperForces(linkLength>meanLinkLength) = linkDampingCoeff*linkLengthDeriv(linkLength>meanLinkLength);

% Total force in each link
totalLinkForceVec = repmat(springForces+damperForces,[2 1]).*linkUnitVecs;

% Force on each node
nodeForceVecs = [totalLinkForceVec(:,1) diff(totalLinkForceVec,1,2) -totalLinkForceVec(:,end)]...
    + vertForce + nodeDragForce;

end

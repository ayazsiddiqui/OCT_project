clear
clc
format compact

%% test signals
dt = 1/100;

amp = 5;
omega = 10;
sim_time = 10;

z_pos = 100;

ini_length = sqrt(z_pos^2 + 2*amp^2);
Vflow = [1;0];

tVec = 0:dt/10:sim_time;

Rn_o_test = timeseries();
Vn_o_test = timeseries();

Rn_o_test.Time = tVec;
Vn_o_test.Time = tVec;

for ii = 1:1
    Rn_o_test.Data(:,ii,:) = [amp*cos(omega*tVec)+ 5*(ii-1);
        z_pos*ones(size(tVec))];
    
    Vn_o_test.Data(:,ii,:) = [-amp*omega*sin(omega*tVec);
        0*ones(size(tVec))];
    
end

%% params
numNode = 4;
tetLength = 95;
tetDia = 0.005;
tetDensity = 1300;
tetYoungs = 4e9;
tetZeta = 0.05;
mT = 2000;
rhoF = 1000;
gravAcc = 9.81;
Cd = 0.5;

R1_o = [0;0];
Rn_o = [0;100];

temp= [linspace(R1_o(1),Rn_o(1),numNode);...
    linspace(R1_o(2),Rn_o(2),numNode)];

temp(:,1) = [];
temp(:,end) = [];

for ii = 1:numNode-2
    init_Ri_o(:,ii) = temp(:,ii);
end

nodePos = [R1_o init_Ri_o Rn_o];
nodeVel = zeros(size(nodePos));

[sumF,nodeMass] = calcTetForce(Vflow,nodePos,nodeVel,tetLength,numNode,...
    tetDia,tetDensity,Cd,tetYoungs,tetZeta,mT,rhoF,gravAcc)

[nodeForceVecs,nodeMass] = nodeForce(Vflow,nodePos,nodeVel,...
    tetLength,mT,tetZeta,tetDia,tetYoungs,Cd,rhoF,tetDensity,gravAcc)

% sim('planarTether_th');
% sim('dummy');



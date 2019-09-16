clear
clc
format compact

% simParam
maxStep = 0.05;

% turbine parameters
Cp = 0.6;
td = 0.15;
Cdamp = 1.5;
Izz = 10;

% state control
A = [0 1;0 0];
B = [0; 1/Izz];
C = [1 0];
D = 0;

Aaug = [A zeros(2,1);C 0];
Baug = [B;0];
H = [0; 0; -1];

eigDes = [-2; -1; -2.5];

K = place(Aaug,Baug,eigDes);

% fluid 
rhoF = 1;

% wake def parameters
nth = 10;
nd = 10;
ky = 0.002;
kz = ky;
CT = 0.9;

upstreamTurbPos = [0;0;0];
downstreamTurbPos = td*[7;0;0];
upstreamTurbFlow = 4;
upstreamTurbYawSpRange = [-30 30]*pi/180;
rngSeed = 1;

% GP sampling time interval
gpSampleInt = 10;
numSample = 10;
simTimeGp = gpSampleInt*numSample;

sim('trainGP_th')

% store training data
trainDsgn = reshape(trainDsgns.Data,1,[],1);
trainFval = reshape(trainFval.Data,[],1,1);

% bayesian ascent parameters
noiseVar = 0.05;
designLimits = upstreamTurbYawSpRange;
gamma = 0.01;
beta = 1.1;
iniTauPerc = 0.05;

simTimeBa = 400;
dt = 20;
maxIter = ceil(simTimeBa/dt);

sim('bayesian_th')







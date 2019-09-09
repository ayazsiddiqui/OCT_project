clear
clc
format compact

% turbine parameters
Cp = 0.6;
td = 0.15;
Cdamp = 1.5;
MI = 10;

% controller parameters
Kp = 0.5;
Ki = 0;
Kd = 4;
tau = 0.5;

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
noiseVar = 0.005;
designLimits = upstreamTurbYawSpRange;
gamma = 0.01;
beta = 1.1;
iniTauPerc = 0.05;

simTimeBa = 600;
dt = 20;
maxIter = ceil(simTimeBa/dt);

sim('bayesian_th')







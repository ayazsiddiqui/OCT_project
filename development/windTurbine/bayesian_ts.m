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
ky = 0.003;
kz = ky;
CT = 0.9;

upstreamTurbPos = [0;0;0];
downstreamTurbPos = td*[7;0;0];
upstreamTurbFlow = 4;
upstreamTurbYawSpRange = [0 50]*pi/180;
rngSeed = 1;

% GP sampling time interval
gpSampleInt = 10;
numSample = 40;
simTimeGp = gpSampleInt*numSample;
trainPts = linspace(upstreamTurbYawSpRange(1),upstreamTurbYawSpRange(2),numSample);
tval = linspace(0,simTimeGp,numSample);
trainSPs = timeseries(trainPts,tval);


sim('trainGP_th')

% store training data
trainDsgn = reshape(trainDsgns.Data,1,[],1);
trainFval = reshape(trainFval.Data,[],1,1);

% bayesian ascent parameters
gp = gaussianProcess(1,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.0005;


designLimits = upstreamTurbYawSpRange;
gamma = 0.01;
beta = 1.1;
iniTauPerc = 0.05;

simTimeBa = 400;
dt = 20;
maxIter = ceil(simTimeBa/dt);

sim('bayesian_th');

%%
time = totPow.time;
Pt = squeeze(totPow.Data);

f1 = figure(1);
plot(time,Pt,'linewidth',1.2)
grid on
hold on
xlabel('Time (s)')
ylabel('$P/P_{max}$')
title('Normalized total power')
saveas(f1,'tPow.png');

f2 = figure(2);
yawC = squeeze(yawCom.Data);
plot(time,yawC*180/pi,'k--','linewidth',1.2)
hold on
grid on
plot(time,squeeze(yawVal.Data)*180/pi,'linewidth',1.2)
xlabel('Time (s)')
ylabel('Angle (deg)')
title('Yaw angle')
legend('$\psi_{SP}$','$\psi$')
saveas(f2,'yawVal.png')








% clear
clear except
clc
format compact
close all

% rngSeed = randi([0,100],1);
rngSeed = 10;
rng(rngSeed);

%% test class
gp = gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';
gp.acquisitionFunctionName = 'expectedImprovement';
gp = gp.buildKernel;
gp = gp.buildAcquitisionFn;

if strcmpi(gp.acquisitionFunctionName,'upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 2;
end

nSamp = 50;
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];
trainDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

trainFval = gp.objectiveFunction(trainDsgns);
noiseVar = 0.005;
gp.kernel.noiseVariance = noiseVar;

%% train GP
% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

%% formulate bayesian ascent
iniTau = 0.05*ones(gp.noInputs,1)*(xMax-xMin);
gamma = 0.01;
beta = 1.1;

maxIter = 10;

iniPt = ((xMax-xMin).*rand(2,1) + xMin);

sim('test')

% trainDsgns,trainFval,trainOpHyp,noiseVar,designLimits
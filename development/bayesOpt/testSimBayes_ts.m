
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
gp.kernel.noiseVariance = 0.005;

%% train GP
% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

sim('test')


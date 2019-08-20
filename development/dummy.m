% dummy

% clear
clear all
clc
format compact
close all

rng('default');
rng(1);

%% test class
gp = OPT.gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';
gp.acquisitionFunction = 'upperConfidenceBound';

nSamp = 20;
xMin = -5; xMax = 5;
testDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

trainFval = gp.objectiveFunction(testDsgns);
gp.getkernel;


% optimize hyper parameters
initialGuess = 1*rand(2+gp.noInputs,1);
opHyp = gp.optimizeHyperParameters(testDsgns,initialGuess);

gp.kernel.covarianceAmp = opHyp(1);
gp.kernel.noiseVariance = opHyp(2);
gp.kernel.lengthScale = opHyp(3:end);

tstCovMat = gp.buildCovarianceMatrix(testDsgns,testDsgns);

%% posterior
% posterior
nPost = 100;
postDsgns = ((xMax-xMin).*rand(2,nPost) + xMin);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,testDsgns,tstCovMat,trainFval);


%% plot results
figure(1)
scatter3(postDsgns(1,:),postDsgns(2,:),gp.objectiveFunction(postDsgns),'k');
hold on
scatter3(testDsgns(1,:),testDsgns(2,:),gp.objectiveFunction(testDsgns)','+b');

% estimations
scatter3(postDsgns(1,:),postDsgns(2,:),predMean,'r+');

grid on
% scatter3(postDsgns,upperLimitMean,'-.m')
% scatter3(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean')

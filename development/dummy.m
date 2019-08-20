% dummy

% clear
clear all
clc
format compact
close all

rng('default');
rng(4);

%% test class
gp = OPT.gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';

nSamp = 20;
xMin = -5; xMax = 5;
x1 = ((xMax-xMin).*rand(1,nSamp) + xMin);
x2 = ((xMax-xMin).*rand(1,nSamp) + xMin);
testDsgns = [x1;x2];

trainFval = gp.objectiveFunction(testDsgns);

gp.getkernel

gp.kernel.covarianceAmp = 1;
gp.kernel.noiseVariance = 0;
gp.kernel.lengthScale = [1;1];

Kmat = gp.buildCovarianceMatrix(testDsgns,testDsgns);

logLike = gp.calcLogLikelihood(testDsgns,'covarianceAmp',1,'noiseVariance',0,'lengthScale',[6;1]);

initialGuess = 1*rand(2+gp.noInputs,1);

val = gp.optimizeHyperParameters(testDsgns,initialGuess);

gp.kernel.covarianceAmp = val(1);
gp.kernel.noiseVariance = val(2);
gp.kernel.lengthScale = val(3:end);

%%
nPost = 500;
x1Post = ((xMax-xMin).*rand(1,nPost) + xMin);
x2Post = ((xMax-xMin).*rand(1,nPost) + xMin);
postDsgns = [x1Post;x2Post];

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,testDsgns,trainFval);

%% plot results
figure(1)
scatter3(postDsgns(1,:),postDsgns(2,:),gp.objectiveFunction(postDsgns),'k');
hold on
scatter3(testDsgns(1,:),testDsgns(2,:),gp.objectiveFunction(testDsgns)','+b');

% estimations
scatter3(postDsgns(1,:),postDsgns(2,:),predMean,'r');
grid on
% scatter3(postDsgns,upperLimitMean,'-.m')
% scatter3(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean')

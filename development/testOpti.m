
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

nSamp = 1;
xMin = -5; xMax = 5;
testDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

trainFval = gp.objectiveFunction(testDsgns);
gp.getkernel;

%% formulate bayesian ascent
percDev = 0.1;

for ii = 1:10
    if ii == 1
        testDsgns = testDsgns;
        trainFval = trainFval;
        initialGuess = 1*rand(2+gp.noInputs,1);
        iniPt = testDsgns(:,randi(size(testDsgns,2)));
    else
        testDsgns = [testDsgns, optPt];
        trainFval = gp.objectiveFunction(testDsgns);
        initialGuess = opHyp;
        iniPt = optPt;
    end
    
    % step 1: optimizie hyper parameters
    opHyp = gp.optimizeHyperParameters(testDsgns,initialGuess);
    gp.kernel.covarianceAmp = opHyp(1);
    gp.kernel.noiseVariance = opHyp(2);
    gp.kernel.lengthScale = opHyp(3:end);
    
    % step 2: construct GP model
    tstCovMat = gp.buildCovarianceMatrix(testDsgns,testDsgns);
    
    % select next input
    xLims = [iniPt - percDev*abs(xMin), iniPt + percDev*xMax];
    aboveLim = xLims > xMax;
    xLims(aboveLim) = xMax;
    
    belowLim = xLims < xMin;
    xLims(belowLim) = xMin;
    
    [optPt,EImax] = gp.maximizeAcquisitionFunction(testDsgns,tstCovMat,trainFval,iniPt,xLims,'explorationFactor',500);
    
end

%% final point
[ma,im] = max(trainFval);
finDsgn = testDsgns(:,im);

%% posterior
% posterior
nPost = 25;
postDsgns = ((xMax-xMin).*rand(2,nPost) + xMin);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,testDsgns,tstCovMat,trainFval);

%% plot results
figure(1)
scatter3(postDsgns(1,:),postDsgns(2,:),gp.objectiveFunction(postDsgns),'k');
hold on
scatter3(testDsgns(1,:),testDsgns(2,:),gp.objectiveFunction(testDsgns)','+b');

% estimations
scatter3(postDsgns(1,:),postDsgns(2,:),predMean,'r+');

% plot final point
scatter3(finDsgn(1,:),finDsgn(2,:),ma,'m*');

grid on
% scatter3(postDsgns,upperLimitMean,'-.m')
% scatter3(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean','Final design')
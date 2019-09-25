% dummy

% clear
clear all
clc
format compact
close all

rng('default');
rng(1);

%% test class
gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
gp.kernel.noiseVariance = 0.001;

nSamp = 50;
xMin = -5; xMax = 5;
testDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);


%% objective functions
% % % Park example 1
% objF = @(X)-((X(1,:).^2 + X(2,:).^2)./50) + 1;
% % % Park example 2
% objF = @(X) 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
%     +0.75*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
% % https://www.hindawi.com/journals/mpe/2013/948303/ example
objF = @(X) exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
    exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
    2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
    1.5.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));

trainFval = objF(testDsgns)';

% optimize hyper parameters
initialGuess = 1*rand(1+gp.noInputs,1);
opHyp = gp.optimizeHyperParameters(testDsgns,trainFval,initialGuess);

gp.kernel.covarianceAmp = opHyp(1);
gp.kernel.lengthScale = opHyp(2:end);

tstCovMat = gp.buildCovarianceMatrix(testDsgns,testDsgns);

%% posterior
% posterior
nPost = 500;
postDsgns = ((xMax-xMin).*rand(2,nPost) + xMin);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,testDsgns,tstCovMat,trainFval);


%% plot results
figure(1)
scatter3(postDsgns(1,:),postDsgns(2,:),objF(postDsgns),'k');
hold on
scatter3(testDsgns(1,:),testDsgns(2,:),objF(testDsgns)','+b');

% estimations
scatter3(postDsgns(1,:),postDsgns(2,:),predMean,'r+');

grid on
% scatter3(postDsgns,upperLimitMean,'-.m')
% scatter3(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean')

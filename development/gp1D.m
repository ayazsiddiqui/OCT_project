% dummy

% clear
clear 
clc
format compact
close all

rng('default');
rng(1);

%% test class
gp = gaussianProcess(1,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
gp.kernel.noiseVariance = 5e-5;

nSamp = 8;
xMin = 0; xMax = 5;
trainDsgns = ((xMax-xMin).*rand(1,nSamp) + xMin);


%% objective functions
% % % Park example 1
% objF = @(X)-((X(1,:).^2 + X(2,:).^2)./50) + 1;
% % % Park example 2
% objF = @(X) 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
%     +0.75*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
% % https://www.hindawi.com/journals/mpe/2013/948303/ example
% objF = @(X) exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
%     exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
%     2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
%     1.5.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));

objF = @(x) 0.25*x.*sin(2*x).*exp(-0.5*x) - 0.125*x.^2.*exp(-x) + 0.12;

trainFval = objF(trainDsgns)';

% optimize hyper parameters
initialGuess = 1*rand(1+gp.noInputs,1);
opHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);
% opHyp = [0.025;0.25];

gp.kernel.covarianceAmp = opHyp(1);
gp.kernel.lengthScale = opHyp(2:end);

trainCovMat = gp.buildCovarianceMatrix(trainDsgns,trainDsgns);

%% posterior
% posterior
nPost = 500;
postDsgns = linspace(xMin,xMax,nPost);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,trainDsgns,trainCovMat,trainFval);
UCB = predMean + 2*predVar;
LCB = predMean - 2*predVar;


%% plot results
f1 = figure(1);
plot(postDsgns(1,:),objF(postDsgns),'k','linewidth',0.9);
hold on
xlabel('$x$')
ylabel('$f(x)$')
plot(trainDsgns(1,:),trainFval,'ob','linewidth',1);

% estimations
plot(postDsgns(1,:),predMean,'r-','linewidth',1);
plot(postDsgns(1,:),UCB,'m--','linewidth',0.9);
plot(postDsgns(1,:),LCB,'m--','linewidth',0.9);

grid on
% scatter3(postDsgns,upperLimitMean,'-.m')
% scatter3(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean','Confidence bounds','Location','best')
title(['$\sigma_{0}$ = ',sprintf('%0.3f', opHyp(1)),', $\theta$ = ',sprintf('%0.3f',opHyp(2)),...
    ', $\sigma_{n}$ = ',sprintf('%.0e', gp.kernel.noiseVariance)]);

saveas(f1,'testGP1D.png');


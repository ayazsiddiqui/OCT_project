
% clear
clear except
clc
format compact
close all

% rngSeed = randi([0,100],1);
rngSeed = 46;
rng(rngSeed);

%% test class
% gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','expectedImprovement');

gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end


nSamp = 40;
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];
trainDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

gp.kernel.noiseVariance = 1*0.05;

%% objective functions
% % % Park example 1
objF = @(X)-((X(1,:).^2 + X(2,:).^2)./50) + 1;
% % % Park example 2
% objF = @(X) 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
%     +0.75*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
% % https://www.hindawi.com/journals/mpe/2013/948303/ example
% objF = @(X) exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
%     exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
%     2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
%     1.5.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));

trainFval = objF(trainDsgns)';

%% train GP
% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

%% formulate bayesian ascent
iniTauPerc = 0.2;
iniTau = iniTauPerc*ones(gp.noInputs,1)*(xMax-xMin);
gamma = 0.01;
beta = 1.1;

iniPt = ((xMax-xMin).*rand(2,1) + xMin);

simTime = 10;
dt = 2;
maxIter = ceil(simTime/dt);
sim('test')

% [sol,gp] = gp.bayesianAscent(trainDsgns,trainFval,trainOpHyp,iniPt,designLimits,iniTau,gamma,beta,maxIter);
% 
% %% final point
% finPts = sol.finPts;
% finDsgn = finPts(:,end);
% testOpHyp = sol.testOpHyp;
% finFval = sol.finFval';



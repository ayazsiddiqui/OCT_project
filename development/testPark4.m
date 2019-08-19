% a proper test based on parks paper

% clear
clear all
clc
format compact
% close all
% repeated random numbers
rng('default');
rng(4);

%% start
% objective function
objF = @(X) -((X(1,:).^2 + X(1,:).^2)./50) + 1;

nSamp = 1;
xMin = -5; xMax = 5;
x1 = ((xMax-xMin).*rand(1,nSamp) + xMin);
x2 = ((xMax-xMin).*rand(1,nSamp) + xMin);
testDsgns = [x1;x2];
testY = objF(testDsgns)';

%% optimize hyper parameters
ini_sigma0 = 0.05;
ini_sigmaE = 0.025;
% initial theta
a = 0;
b = 1;
ini_theta = ((b-a).*rand(size(testDsgns,1),1) + a);
ini_hyperParam = [ini_sigma0;ini_sigmaE;ini_theta];
% constraints on optimization
A = [];b = [];Aeq = [];beq = [];
lb = 1e-2*ones(2+numel(ini_theta),1);
ub = 100*ones(2+numel(ini_theta),1);

% optimize using fmincon
% [optHyper,fval] = fmincon(@(hyperParam) calcLogLogLikelihood(testY,testDsgns,hyperParam),ini_hyperParam,A,b,Aeq,beq,lb,ub);

%%
% optHyper(1) = 1*0.125;
% optHyper(2) = 20;
% optHyper(3) = 2*0.18;
optHyper = [
    0.6584
    0.0100
    5.2535
   99.9998];

% post designs
nPost = 2;
x1Post = ((xMax-xMin).*rand(1,nPost) + xMin);
x2Post = ((xMax-xMin).*rand(1,nPost) + xMin);
postDsgns = [x1Post;x2Post];

[muD,Var] = posteriorCalc(testDsgns,testY,optHyper,postDsgns);

upperLimitMean = muD + 2*Var;
lowerLimitMean = muD - 2*Var;

x1Fin = ((xMax-xMin).*rand(1,1) + xMin);
x2Fin = ((xMax-xMin).*rand(1,1) + xMin);
XFin = [x1Fin;x2Fin];

EI = calcExpectedImprovement(testDsgns,testY,optHyper,XFin);


%%
A = [];b = [];Aeq = [];beq = [];
lb = -5*ones(2,1);
ub = 5*ones(2,1);

[XOpt,fval] = fmincon(@(XFin) calcExpectedImprovement(testDsgns,testY,optHyper,XFin),XFin,A,b,Aeq,beq,lb,ub);


%% plot results
figure(1)
scatter3(postDsgns(1,:),postDsgns(2,:),objF(postDsgns),'k');
hold on
scatter3(testDsgns(1,:),testDsgns(2,:),objF(testDsgns)','+b');

% estimations
scatter3(postDsgns(1,:),postDsgns(2,:),muD,'r');
grid on
% scatter3(postDsgns,upperLimitMean,'-.m')
% scatter3(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean')






% a proper test based on parks paper

% clear
clear all
clc
format compact
close all
% repeated random numbers
rng('default');
rng(1);

%% start
% objective function
objF = @(X) 0.25*X(1,:).*sin(2*X(1,:)).*exp(-0.5*X(1,:)) - 0.125*X(1,:).^2.*exp(-X(1,:)) + 0.12;

nSamp = 20;
% testDsgns = ((5-0).*rand(nSamp,1) + 0)';
testDsgns = linspace(0,5,nSamp);
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
lb = zeros(2+numel(ini_theta),1);
ub = 100*ones(2+numel(ini_theta),1);

% optimize using fmincon
[optHyper,fval] = fmincon(@(hyperParam) calcLogLogLikelihood(testY,testDsgns,hyperParam),ini_hyperParam,A,b,Aeq,beq,lb,ub);

%%
% optHyper(1) = 5*0.125;
% optHyper(2) = 0.009;
% optHyper(3) = 3*0.18;

% post designs
nPost = 50;
postDsgns = linspace(0,5,nPost);

CovMatVecTranspose = buildCovMat(testDsgns,postDsgns,'covAmplitude',optHyper(1),'noiseVariance',optHyper(2),'lengthScale',optHyper(3:end));
testCovMat = buildCovMat(testDsgns,testDsgns,'covAmplitude',optHyper(1),'noiseVariance',optHyper(2),'lengthScale',optHyper(3:end));

% CovMatVecTranspose
% muD = (CovMatVecTranspose'/testCovMat)*testY';
% Var = buildCovMat(postDsgns,postDsgns,'covAmplitude',optHyper(1),'noiseVariance',optHyper(2),'lengthScale',optHyper(3:end));

%% mean and variance
muD = NaN(nPost,1);
Var = NaN(nPost,1);
for ii = 1:nPost
    x= 1;    
    muD(ii,1) = (CovMatVecTranspose(:,ii)'/testCovMat)*testY;
    Var(ii,1) = buildCovMat(postDsgns(:,ii),postDsgns(:,ii),'covAmplitude',optHyper(1),'noiseVariance',optHyper(2),'lengthScale',optHyper(3:end)) - (CovMatVecTranspose(:,ii)'/testCovMat)*CovMatVecTranspose(:,ii);
end  

upperLimitMean = muD + 2*Var;
lowerLimitMean = muD - 2*Var;


%% plot results
figure
plot(postDsgns,objF(postDsgns),'--k');
hold on
scatter(testDsgns,objF(testDsgns)','+b');

% estimations
plot(postDsgns,muD,'-r');
grid on
plot(postDsgns,upperLimitMean,'-.m')
plot(postDsgns,lowerLimitMean,'-.m')

legend('objF','Sampled points','Estimated mean','Upper Conf. int.','Lower Conf. int.')






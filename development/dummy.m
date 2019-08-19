% dummy

% clear
clear all
clc
format compact
close all

testDsgns = linspace(0,5,10);
testDsgns = [testDsgns;testDsgns];

covMat = buildCovMat(testDsgns,testDsgns,'covAmplitude',0.2,'noiseVariance',0,'lengthScale',[1;1]);


objF = @(X) (0.25*X(1,:).*sin(2*X(1,:)).*exp(-0.5*X(1,:)) - 0.125*X(1,:).^2.*exp(-X(1,:)) + 0.12)';

x = 1;
% op = calcLogLogLikelihood(objF(testDsgns),testDsgns,[0.2,0,1]);

%% test class
gp = OPT.gaussianProcess;

gp.setNoInputs(2);
gp.setKernel('squaredExponential');

nSamp = 4;
xMin = -5; xMax = 5;
x1 = ((xMax-xMin).*rand(1,nSamp) + xMin);
x2 = ((xMax-xMin).*rand(1,nSamp) + xMin);
testDsgns = [x1;x2];

gp.objectiveFunction(testDsgns)

Kmat = gp.buildCovarianceMatrix(testDsgns,testDsgns,'covAmplitude',0.5,'noiseVariance',0,'lengthScale',[1;1]);





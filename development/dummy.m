% dummy

% clear
clear all
clc
format compact
close all

testDsgns = linspace(0,5,10);
testDsgns = [testDsgns;testDsgns];

covMat = buildCovMat(testDsgns,testDsgns,'covAmplitude',0.2,'noiseVariance',0,'lengthScale',[1;1]);


% objF = @(X) (0.25*X(1,:).*sin(2*X(1,:)).*exp(-0.5*X(1,:)) - 0.125*X(1,:).^2.*exp(-X(1,:)) + 0.12)';

x = 1;
% op = calcLogLogLikelihood(objF(testDsgns),testDsgns,[0.2,0,1]);
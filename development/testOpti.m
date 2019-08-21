
% clear
clear all
clc
format compact
close all

rng('default');
rng(5);

% good runs:
% rng = 1,2,15,6

%% test class
gp = OPT.gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';
gp.acquisitionFunction = 'upperConfidenceBound';

nSamp = 1;
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];
testDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

trainFval = gp.objectiveFunction(testDsgns);
gp.getkernel;

%% formulate bayesian ascent
iniTau = 0.05*ones(gp.noInputs,1)*(xMax-xMin);
beta = 1.05;

for ii = 1:20
    if ii == 1
        testDsgns = testDsgns;
        trainFval = trainFval;
        initialGuess = 1*rand(2+gp.noInputs,1);
        iniPt = testDsgns(:,randi(size(testDsgns,2)));
        tau = iniTau;
    else
        testDsgns = [testDsgns, optPt];
        trainFval = gp.objectiveFunction(testDsgns);
        initialGuess = opHyp;
        iniPt = optPt;
        if trainFval(end) >= (1/ii)*(max(trainFval)-trainFval(1))
            tau = beta*tau;
        else
            tau = iniTau;
            pause(5);
            fprintf('Bounds reset to initial bounds at iteration %d',ii)
        end
        
    end
    
    % step 1: optimizie hyper parameters
    opHyp = gp.optimizeHyperParameters(testDsgns,initialGuess);
    gp.kernel.covarianceAmp = opHyp(1);
    gp.kernel.noiseVariance = opHyp(2);
    gp.kernel.lengthScale = opHyp(3:end);
    
    % step 2: construct GP model
    tstCovMat = gp.buildCovarianceMatrix(testDsgns,testDsgns);
    
    % select next input
    resetSwitch = 0;
    xLims = gp.calDesignBounds(iniPt,tau,designLimits);
    
    [optPt,EImax] = gp.maximizeAcquisitionFunction(testDsgns,tstCovMat,trainFval,iniPt,xLims,'explorationFactor',2);
    
end

%% final point
[ma,im] = max(trainFval);
finDsgn = testDsgns(:,im);

%% posterior
% posterior
nPost = 50;
postDsgns = ((xMax-xMin).*rand(2,nPost) + xMin);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,testDsgns,tstCovMat,trainFval);

% %% plot results
% figure(1)
% scatter3(postDsgns(1,:),postDsgns(2,:),gp.objectiveFunction(postDsgns),'k');
% hold on
% scatter3(testDsgns(1,:),testDsgns(2,:),gp.objectiveFunction(testDsgns)','+b');
%
% % estimations
% scatter3(postDsgns(1,:),postDsgns(2,:),predMean,'r+');
%
% % plot final point
% scatter3(finDsgn(1,:),finDsgn(2,:),ma,'m*');
%
% grid on
% % scatter3(postDsgns,upperLimitMean,'-.m')
% % scatter3(postDsgns,lowerLimitMean,'-.m')
%
% legend('objF','Sampled points','Estimated mean','Final design')

%% plot grid
nGrid = 25;
x1Grid = linspace(xMin,xMax,nGrid);
x2Grid = linspace(xMin,xMax,nGrid);

[X1,X2] = meshgrid(x1Grid,x2Grid);

Z = gp.objectiveFunction([X1(:)';X2(:)']);
Z = reshape(Z,nGrid,nGrid);

% figure
figure(2)
contourf(X1,X2,Z)
colorbar
hold on
plot(testDsgns(1,:),testDsgns(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot(finDsgn(1,:),finDsgn(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize',10)


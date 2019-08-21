
% clear
clear all
clc
format compact
close all

rng('default');
rng(8);

% good runs:
% rng = 1,2,15,6

%% test class
gp = OPT.gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';
gp.acquisitionFunction = 'upperConfidenceBound';

nSamp = 20;
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];
trainDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

trainFval = gp.objectiveFunction(trainDsgns);
gp.getkernel;
gp.kernel.noiseVariance = 0.005;

%% train GP
% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,initialGuess);

% step 2: construct GP model
trainCovMat = gp.buildCovarianceMatrix(trainDsgns,trainDsgns);


%% formulate bayesian ascent
iniTau = 0.05*ones(gp.noInputs,1)*(xMax-xMin);
gamma = 0.05;
beta = 1.1;

for ii = 1:10
    if ii == 1
        testDsgns = trainDsgns;
        testFval = trainFval;
        testOpHyp = trainOpHyp;
        iniPt =  ((xMax-xMin).*rand(2,1) + xMin);
        finPts = iniPt;
        finFval = gp.objectiveFunction(iniPt);
        tau = iniTau;
        
    else
        testDsgns = [testDsgns, optPt];
        testFval = [testFval; optFval];
        testOpHyp = gp.optimizeHyperParameters(testDsgns,testOpHyp);
        iniPt = optPt;
        finPts = [finPts,optPt];
        finFval = [finFval; optFval];
        
        if finFval(end) >= gamma*(1/ii)*(max(trainFval)-finFval(1))
            tau = beta*tau;
            
        else
            tau = iniTau;
            pause(5);
            fprintf('Bounds reset to initial bounds at iteration %d',ii)
        end
        
    end
    
    % step 1: optimizie hyper parameters
    gp.kernel.covarianceAmp = testOpHyp(1);
    gp.kernel.lengthScale = testOpHyp(2:end);
    
    % step 2: construct GP model
    tstCovMat = gp.buildCovarianceMatrix(testDsgns,testDsgns);
    
    % select next input
    xLims = gp.calDesignBounds(iniPt,tau,designLimits);
    
    % maximize acquisition function
    [optPt,AQmax] = gp.maximizeAcquisitionFunction(testDsgns,tstCovMat,testFval,iniPt,xLims,'explorationFactor',2.5);
    optFval = gp.objectiveFunction(optPt);
    
end

%% final point
[ma,im] = max(testFval);
finDsgn = testDsgns(:,im);

%% posterior
% posterior
nPost = 50;
postDsgns = ((xMax-xMin).*rand(2,nPost) + xMin);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,testDsgns,tstCovMat,testFval);

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
plot(finPts(1,:),finPts(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot(finDsgn(1,:),finDsgn(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize',10)


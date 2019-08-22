
% clear
clear except -runNo
clc
format compact
% close all

rngSeed = 60;
rng('default');
rng(rngSeed);

%% test class
gp = OPT.gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';
gp.acquisitionFunction = 'upperConfidenceBound';
gp.acquisitionFunction = 'expectedImprovement';

nSamp = 50;
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
gamma = 0.01;
beta = 1.1;

knownMax = 2;

noIter = 1;
goodNess = 0;
% && goodNess < 0.99
while noIter <= 20 
    if noIter == 1
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
        
        if finFval(end) >= gamma*(1/noIter)*(max(trainFval)-finFval(1))
            tau = beta*tau;
            
        else
            tau = iniTau;
%             pause(5);
            fprintf('Bounds reset to initial bounds at iteration %d\n',noIter)
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
    [optPt,AQmax] = gp.maximizeAcquisitionFunction(testDsgns,tstCovMat,testFval,finFval,iniPt,xLims,'explorationFactor',2.5);
    optFval = gp.objectiveFunction(optPt);
    
    % convergence check
    goodNess = optFval/knownMax;
    noIter = noIter + 1;
    
end

%% final point
% [ma,im] = max(finFval);
% finDsgn = finPts(:,im);
finDsgn = optPt;

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
nGrid = 100;
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
plot([finPts(1,:),finDsgn(1,:)],[finPts(2,:),finDsgn(2,:)],'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot(finDsgn(1,:),finDsgn(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize',10)
xlabel('$x_{1}$')
ylabel('$x_{2}$')
title(sprintf('RNG seed = %d',rngSeed))

% % figure
% if exist('runNo','var')
%     runNo = runNo+1;
% else
%     runNo = 1;
% end
% figure(3)
% set(gcf,'Position',[50 50 3.5*560 2*420])
% subplot(2,4,runNo)
% contourf(X1,X2,Z)
% hold on
% plot([finPts(1,:),finDsgn(1,:)],[finPts(2,:),finDsgn(2,:)],'-rs',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',10)
% plot(finDsgn(1,:),finDsgn(2,:),'-rs',...
%     'LineWidth',2,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','m',...
%     'MarkerSize',10)
% title(sprintf('RNG seed = %d',rngSeed))
% axis equal
% xlabel('$x_{1}$')
% xlabel('$x_{2}$')

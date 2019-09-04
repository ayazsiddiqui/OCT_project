
% clear
clear except
clc
format compact
close all

% rngSeed = randi([0,100],1);
rngSeed = 10;
rng(rngSeed);

%% test class
gp = gaussianProcess;

gp.noInputs = 2;
gp.kernelName = 'squaredExponential';
gp.acquisitionFunctionName = 'expectedImprovement';
gp = gp.buildKernel;
gp = gp.buildAcquitisionFn;

if strcmpi(gp.acquisitionFunctionName,'upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 2;
end

nSamp = 50;
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];
trainDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

trainFval = gp.objectiveFunction(trainDsgns);
gp.kernel.noiseVariance = 0.005;

%% train GP
% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

%% formulate bayesian ascent
iniTau = 0.05*ones(gp.noInputs,1)*(xMax-xMin);
gamma = 0.01;
beta = 1.1;

maxIter = 10;

iniPt = ((xMax-xMin).*rand(2,1) + xMin);

[sol,gp] = gp.bayesianAscent(trainDsgns,trainFval,trainOpHyp,iniPt,designLimits,iniTau,gamma,beta,maxIter);


%% final point
finPts = sol.finPts;
finDsgn = finPts(:,end);

%% posterior
% posterior
nPost = 50;
postDsgns = ((xMax-xMin).*rand(2,nPost) + xMin);

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,sol.testDsgns,sol.testCovMat,sol.testFval);

%% plot grid
nGrid = 100;
x1Grid = linspace(xMin,xMax,nGrid);
x2Grid = linspace(xMin,xMax,nGrid);

[X1,X2] = meshgrid(x1Grid,x2Grid);

Z = gp.objectiveFunction([X1(:)';X2(:)']);
Z = reshape(Z,nGrid,nGrid);

% figure
locs = getFigLocations(2*560,420);
figure(2)
set(gcf,'Position',locs(1,:))
subplot(1,2,1)
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
xlabel('$x_{1}$')
ylabel('$x_{2}$')
title(sprintf('EI, RNG seed = %d',rngSeed))

%% run again
gp.acquisitionFunctionName = 'upperConfidenceBound';
gp = gp.buildAcquitisionFn;

if strcmpi(gp.acquisitionFunctionName,'upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 5;
end
[sol,gp] = gp.bayesianAscent(trainDsgns,trainFval,trainOpHyp,iniPt,designLimits,iniTau,gamma,beta,maxIter);

%% final point
finPts = sol.finPts;
finDsgn = finPts(:,end);

% figure
fg = figure(2);
subplot(1,2,2)
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
xlabel('$x_{1}$')
ylabel('$x_{2}$')
title(sprintf('UCB, RNG seed = %d',rngSeed))

%% saveas
% saveas(fg,sprintf('figNo%d.png',rngSeed));



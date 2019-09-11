
% clear
clear except
clc
format compact
close all

% rngSeed = randi([0,100],1);
rngSeed = 40;
rng(rngSeed);

%% test class
gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','expectedImprovement');
gp.kernel.noiseVariance = 0.005;

if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 2;
end

nSamp = 10;
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];
trainDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

%% objective functions
% % % Park example 1
% objF = @(X)-((X(1,:).^2 + X(2,:).^2)./50) + 1;
% % % Park example 2
objF = @(X) 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
    +0.5*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
% % https://www.hindawi.com/journals/mpe/2013/948303/ example
% objF = @(X) exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
%     exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
%     2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
%     2.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));
            
trainFval = objF(trainDsgns);
trainFval = trainFval(:);

%% train GP
% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

%% formulate bayesian ascent
iniTau = 0.1*ones(gp.noInputs,1)*(xMax-xMin);
gamma = 0.01;
beta = 1.1;

maxIter = 20;

iniPt = ((xMax-xMin).*rand(2,1) + xMin);
finPtsEI = iniPt;
iniFval = objF(iniPt);
finFvalEI = iniFval;
opHypEI = [];
tauEI = [];
predMeanEI = [];
predVarEI = [];

for noIter = 1:maxIter
    
    [sol,gp] = gp.bayesianAscent(trainDsgns,trainFval,finPtsEI,finFvalEI,...
        opHypEI,tauEI,designLimits,iniTau,gamma,beta,noIter);
    
    finPtsEI = [finPtsEI sol.optPt];
    finFvalEI = [finFvalEI;objF(sol.optPt)];
    opHypEI = sol.testOpHyp;
    tauEI = sol.tau;
    predMeanEI = [predMeanEI sol.predMean];
    predVarEI = [predVarEI sol.predVar];
    
end

postDsgns = ((xMax-xMin).*rand(2,4) + xMin);

val = gp.calcAcquisitionFunction(postDsgns,max(finFvalEI),sol.testDsgns,sol.testCovMat,sol.testFval);

%% plot grid
nGrid = 100;
x1Grid = linspace(xMin,xMax,nGrid);
x2Grid = linspace(xMin,xMax,nGrid);

[X1,X2] = meshgrid(x1Grid,x2Grid);

Z = objF([X1(:)';X2(:)']);
Z = reshape(Z,nGrid,nGrid);

% figure
fidWid = 340;
locs = getFigLocations(2*(4/3)*fidWid,fidWid);
figure(1)
set(gcf,'Position',locs(1,:))
subplot(1,2,1)
contourf(X1,X2,Z)
colorbar
hold on
plot(finPtsEI(1,:),finPtsEI(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot(finPtsEI(1,end),finPtsEI(2,end),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize',10)
xlabel('$x_{1}$')
ylabel('$x_{2}$')
title(sprintf('EI, RNG seed = %d',rngSeed))

figure(2)
set(gcf,'Position',locs(2,:))
subplot(1,2,1)
plot(1:maxIter,predMeanEI,'-o');
grid on
hold on
xlabel('Iteration number')
ylabel('Predicted mean')
title(sprintf('EI, RNG seed = %d',rngSeed))

figure(3)
set(gcf,'Position',locs(3,:))
subplot(1,2,1)
plot(1:maxIter,predVarEI,'-o');
grid on
hold on
xlabel('Iteration number')
ylabel('Predicted variance')
title(sprintf('EI, RNG seed = %d',rngSeed))

%% run again
%% test class
gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
gp.kernel.noiseVariance = 0.005;

if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 2;
end

%% bayesian ascent with UCB
finPtsUCB = iniPt;
finFvalUCB = iniFval;
opHypUCB = [];
tauUCB = [];
predMeanUCB = [];
predVarUCB = [];


for noIter = 1:maxIter
    
    [sol,gp] = gp.bayesianAscent(trainDsgns,trainFval,finPtsUCB,finFvalUCB,...
        opHypUCB,tauUCB,designLimits,iniTau,gamma,beta,noIter);
    
    finPtsUCB = [finPtsUCB sol.optPt];
    finFvalUCB = [finFvalUCB;objF(sol.optPt)];
    opHypUCB = sol.testOpHyp;
    tauUCB = sol.tau;
    predMeanUCB = [predMeanUCB sol.predMean];
    predVarUCB = [predVarUCB sol.predVar];

end


val = gp.calcAcquisitionFunction(postDsgns,max(finFvalUCB),sol.testDsgns,sol.testCovMat,sol.testFval);

%% figure
figure(1)
subplot(1,2,2)
contourf(X1,X2,Z)
colorbar
hold on
plot(finPtsUCB(1,:),finPtsUCB(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot(finPtsUCB(1,end),finPtsUCB(2,end),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize',10)
xlabel('$x_{1}$')
ylabel('$x_{2}$')
title(sprintf('UCB, RNG seed = %d',rngSeed))

figure(2)
subplot(1,2,2)
plot(1:maxIter,predMeanUCB,'-o');
grid on
hold on
xlabel('Iteration number')
ylabel('Predicted mean')
title(sprintf('UCB, RNG seed = %d',rngSeed))

figure(3)
subplot(1,2,2)
plot(1:maxIter,predVarUCB,'-o');
grid on
hold on
xlabel('Iteration number')
ylabel('Predicted variance')
title(sprintf('UCB, RNG seed = %d',rngSeed))

%% saveas
% saveas(fg,sprintf('figNo%d.png',rngSeed));




% clear
clear
clc
format compact
close all

% rngSeed = randi([0,100],1);
rngSeed = 9;
rng(rngSeed);

%% test class
gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','expectedImprovement');
gp.kernel.noiseVariance = 1*0.05;

nSamp = 12;
xMin = -5; xMax = 5;
x = linspace(xMin,xMax,nSamp);
[x1s,x2s] = meshgrid(x,x);
designLimits = [xMin*[1;1],xMax*[1;1]];
trainDsgns = [x1s(:)';x2s(:)'];
% trainDsgns = ((xMax-xMin).*rand(2,nSamp) + xMin);

%% objective functions
% % % Park example 1
% objF = @(X)-((X(1,:).^2 + X(2,:).^2)./50) + 1;
% % % Park example 2
objF = @(X) 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
    +0.75*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
% % https://www.hindawi.com/journals/mpe/2013/948303/ example
% objF = @(X) exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
%     exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
%     2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
%     1.5.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));


%% train GP
trainFval = objF(trainDsgns);
trainFval = trainFval(:);

%% formulate bayesian ascent
iniTau = 0.2*ones(gp.noInputs,1)*(xMax-xMin);
gamma = 0.01;
beta = 1.1;

iniPt = ((xMax-xMin).*rand(2,1) + xMin);
% iniPt = [5;5];
finPtsEI = iniPt;
iniFval = objF(iniPt);
finFvalEI = iniFval;
opHypEI = [];
tauEI = [];
predMeanEI = [];
predVarEI = [];
AqFnEI = [];
expFac = 1;
maxIter = 5;
predHorizon = 5;
ctrlHorizon = 1;
% 
% [optDsgn,maxF] = particleSwarmOpt(@(x)objF(x),iniPt,designLimits(:,1),designLimits(:,2),...
%     'swarmSize',25,'cognitiveLR',0.4,'socialLR',0.2,'maxIter',20);


for noIter = 1:maxIter
    
    [sol,gp] = gp.mpcBayesianAscent(trainDsgns,trainFval,finPtsEI,finFvalEI,...
        opHypEI,tauEI,designLimits,iniTau,gamma,beta,noIter,predHorizon,ctrlHorizon);
    
    finPtsEI = [finPtsEI sol.optPt];
    finFvalEI = [finFvalEI;objF(sol.optPt)'];
    opHypEI = sol.testOpHyp;
    tauEI = sol.tau;
    predMeanEI = [predMeanEI sol.mpcPredMean];
    predVarEI = [predVarEI sol.mpcPredVar];
    AqFnEI = [AqFnEI sol.optAq];
    
end

postDsgns = ((xMax-xMin).*rand(2,500) + xMin);
[postPredMean,PredVar] = gp.calcPredictiveMeanAndVariance(postDsgns,sol.testDsgns,sol.testCovMat,sol.testFval);



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
f1 = figure(1);
set(gcf,'Position',locs(1,:))
subplot(1,2,1)
contourf(X1,X2,Z)
colormap jet
colorbar
hold on
plot(finPtsEI(1,:),finPtsEI(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)
plot(finPtsEI(1,1),finPtsEI(2,1),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',10)
plot(finPtsEI(1,end),finPtsEI(2,end),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','m',...
    'MarkerSize',10)
xlabel('$x_{1}$')
ylabel('$x_{2}$')
title(sprintf('EI, RNG seed = %d',rngSeed))

f2 = figure(2);
set(gcf,'Position',locs(2,:))
subplot(1,2,1)
for ii = 1:maxIter
    plot(1:predHorizon,predMeanEI(:,ii),'-o');
    hold on
end
grid on
hold on
xlabel('Predcition horizon')
ylabel('Predicted mean')
title(sprintf('EI, RNG seed = %d',rngSeed))

f3 = figure(3);
set(gcf,'Position',locs(3,:))
subplot(1,2,1)
for ii = 1:maxIter
    plot(1:predHorizon,predVarEI(:,ii),'-o');
    hold on
end
grid on
xlabel('Predcition horizon')
ylabel('Predicted variance')
title(sprintf('EI, RNG seed = %d',rngSeed))

figure(4);
set(gcf,'Position',locs(4,:))
subplot(1,2,1)
for ii = 1:maxIter
    plot(ii*ones(1,predHorizon),AqFnEI(:,ii),'-o');
    hold on
end
grid on
xlabel('Iteration number')
ylabel('Objective functional')
title(sprintf('EI, RNG seed = %d',rngSeed))

%% run again
%% test class
gp = gaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
gp.kernel.noiseVariance = 0.005;

if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = expFac;
end

%% bayesian ascent with UCB
finPtsUCB = iniPt;
finFvalUCB = iniFval;
opHypUCB = [];
tauUCB = [];
predMeanUCB = [];
predVarUCB = [];
AqFnUCB = [];

for noIter = 1:maxIter
    
    [sol,gp] = gp.mpcBayesianAscent(trainDsgns,trainFval,finPtsEI,finFvalEI,...
        opHypEI,tauEI,designLimits,iniTau,gamma,beta,noIter,predHorizon,ctrlHorizon);
    
    finPtsUCB = [finPtsUCB sol.optPt];
    finFvalUCB = [finFvalUCB;objF(sol.optPt)'];
    opHypUCB = sol.testOpHyp;
    tauUCB = sol.tau;
    AqFnUCB = [AqFnUCB  sol.optAq];
    predMeanUCB = [predMeanUCB sol.mpcPredMean];
    predVarUCB = [predVarUCB sol.mpcPredVar];
    
end

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
plot(finPtsUCB(1,1),finPtsUCB(2,1),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
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
for ii = 1:maxIter
    plot(1:predHorizon,predMeanUCB(:,ii),'-o');
    hold on
end
grid on
hold on
xlabel('Predcition horizon')
ylabel('Predicted mean')
title(sprintf('UCB, RNG seed = %d',rngSeed))

figure(3)
subplot(1,2,2)
for ii = 1:maxIter
    plot(1:predHorizon,predVarUCB(:,ii),'-o');
    hold on
end
grid on
hold on
xlabel('Predcition horizon')
ylabel('Predicted variance')
title(sprintf('UCB, RNG seed = %d',rngSeed))

figure(4)
set(gcf,'Position',locs(4,:))
subplot(1,2,2)
for ii = 1:maxIter
    plot(ii*ones(1,predHorizon),AqFnUCB(:,ii),'-o');
    hold on
end
grid on
xlabel('Iteration number')
ylabel('Objective functional')
title(sprintf('UCB, RNG seed = %d',rngSeed))

figure(5)
set(gcf,'Position',locs(5,:))
hold on
view(-30,45)
surf(X1,X2,Z)
hold on
scatter3(postDsgns(1,:),postDsgns(2,:),postPredMean)
xlabel('$x_{1}$')
ylabel('$x_{2}$')
zlabel('$ObjF$')



%% saveas
saveas(f1,sprintf('contour%d.png',rngSeed));
saveas(f2,sprintf('predMean%d.png',rngSeed));
saveas(f3,sprintf('predVAr%d.png',rngSeed));



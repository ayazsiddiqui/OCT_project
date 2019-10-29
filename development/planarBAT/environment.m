clear
clc
format compact
rng(2);

%% environment
hMax = 1000;
hMin = 100;
hStep = 100;
heights = hMin:hStep:hMax;
heights = heights(:);
meanFlow = 10;

tVec = 0:5:60;
timeInSec = 60*tVec;

stdDev = 0.8;
timeScale = 30;
heightScale = 200;
windSpeedOut = genWindv2(heights,heightScale,tVec,timeScale,stdDev);
[hMesh,tMesh] = meshgrid(heights,tVec);

nS = length(tVec);


for ii = 1:nS
    Flows(:,:,ii) = meanFlow*(1 + windSpeedOut(:,ii));
        
end

objF = @(h,t,Flows,hMesh,tMesh) interp2(hMesh,60*tMesh,reshape(Flows,size(hMesh)),h,t);

% plot

for ii = 1:nS
    if ii == 1
        grid on
        hold on
        xlabel('Flow speed (m/s)')
        ylabel('Altitude (m)')
        xlim(max(ceil(abs(Flows)),[],'all')*[0 1]);
    else
        delete(pF)
    end
    
    pF = plot(Flows(:,:,ii),heights,'k');
    title(sprintf('Time = %0.1f min',tVec(ii)));
%     waitforbuttonpress
end
        

%% train GP
gp = timeDepGaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.05;

tVals = reshape(transpose(timeInSec'.*ones(length(timeInSec),length(heights))),[],1)';
trainDsgns = [repmat(heights',1,nS);tVals];
trainFval = Flows(:);

% step 1: optimize hyper parameters
initialGuess = rand(1+gp.noInputs,1);
trainOpHyp = gp.optimizeHyperParameters(trainDsgns,trainFval,initialGuess);

gp.kernel.covarianceAmp = trainOpHyp(1);
gp.kernel.lengthScale = trainOpHyp(2:end);

trainCovMat = gp.buildCovarianceMatrix(trainDsgns,trainDsgns);

%% formulate bayesian ascent
iniTau = 0.1*ones(gp.noInputs-1,1)*(hMax-hMin);
gamma = 0.01;
beta = 1.1;
designLimits = [hMin hMax];

iniPt = [500;tVals(end)];
finPtsEI = iniPt;
iniFval = objF(iniPt(1),iniPt(2),windSpeedOut,hMesh,tMesh);
finFvalEI = iniFval;
opHypEI = [];
tauEI = [];
predMeanEI = [];
predVarEI = [];
AqFnEI = [];
maxIter = 5;
predSteps = 5;
timeStep = 100;


tic
for noIter = 1:maxIter
    
    [sol,gp] = gp.mpcBayesianAscent(trainDsgns,trainFval,finPtsEI,finFvalEI,...
        opHypEI,tauEI,designLimits,iniTau,gamma,beta,noIter,predSteps,timeStep);
    
    finPtsEI = [finPtsEI sol.optPt];
    finFvalEI = [finFvalEI;polyval(pt,sol.optPt(1))'];
    opHypEI = [opHypEI sol.testOpHyp];
    tauEI = sol.tau;
    predMeanEI = [predMeanEI sol.mpcPredMean];
    predVarEI = [predVarEI sol.mpcPredVar];
    AqFnEI = [AqFnEI sol.optAq];
    
end

toc


%% posterior
% posterior
nPost = 10;
postAlts = ((heights(end)-heights(1)).*rand(1,nPost) + heights(1));
postAlts = reshape(heights,1,[]);
postTime = 2500*ones(size(postAlts));

postDsgns = [postAlts;postTime];

[predMean,predVar] = gp.calcPredictiveMeanAndVariance(postDsgns,trainDsgns,trainCovMat,trainFval);


%% plot
figure(1)
xlabel('Flow speed (m/s)')
ylabel('Height (m)')
hold on
grid on

for ii = 1:nS
    if ii>1
        delete(p1)
    end
    p1 = plot(Flows(:,:,ii),heights);
    title(sprintf('Time = %0.1f s',FlowInt(ii)))
%     waitforbuttonpress
end





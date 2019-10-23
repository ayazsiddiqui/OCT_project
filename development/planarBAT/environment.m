clear
clc
format compact
rng(2);

%% environment
hMax = 1000;
hMin = 0;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;

% fit a polynomial to dummy data
ht = hMin:200:hMax;
ht = ht(:);
ft2 = [5;5.5;6.5;7;6;5.5];

pt = polyfit(ht,ft2,3);

Flows = polyval(pt,heights);

% simTime in minutes
simTime = 30;
FlowInt = 60*(0:5:simTime);
nS = length(FlowInt);

for ii = 1:nS
    
    for jj = 1:length(pt)
        rdNum = rand;
        if rdNum < 0.33
            mut = -1;
        elseif rdNum >= 0.33 && rdNum < 0.66
            mut = 0;
        else
            mut = 1;
        end
        pt(jj) = pt(jj)*(1 + mut*0.1);
    end
    Flows(:,:,ii) = polyval(pt,heights);
end

%% train GP
gp = timeDepGaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.05;

tVals = reshape(transpose(FlowInt'.*ones(length(FlowInt),length(heights))),[],1)';
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
iniFval = polyval(pt,iniPt(1));
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





clear
clc
format compact

%% initialization
% construct an instance of the RGP class
rgp = RGP(1);

% set grid size
size_grid = 201;

% make basis vector space
basisVecs = linspace(0,20,size_grid)*(pi/180);
% populate initial space with dummy values for the latent function
dummyVals = rand(size_grid,1);

% know values
desPtSeq = [0 0 1 2 3 4 15 5 7 6 8 12 8 8 8 8 8]*(pi/180);
PerfmIndxLib = [0.3 0.28 0.22 0.35 0.18 0.4 0.33 0.42 0.24 0.26 0.19 0.25 0.24 0.29 0.3 0.31 0.19]';

% set values of hyper parameters
hyp.funcVar = 1;         % variance of latent function
hyp.lengthScale = 0.01;  % length scale
hyp.noiseVar = 0.12 ;    % noise variance
hypParams = [hyp.funcVar;hyp.lengthScale;hyp.noiseVar];

% test function that calculates covariance and mean
[K,M] = rgp.buildCovMatAndMeanVec(basisVecs,hypParams);

% test RGP
initVals.basisVec = basisVecs;
[initVals.initCovMat,initVals.initMeanVec] = ...
    rgp.buildCovMatAndMeanVec(basisVecs,hypParams);
initVals.initInvCovMat = inv(initVals.initCovMat + (1e-5)*eye(size_grid));

cGt_1 = initVals.initCovMat;
muGt_1 = initVals.initMeanVec;

for stepNo = 1:length(desPtSeq)
    [muGt,cGt] = rgp.rgpRegression(desPtSeq(stepNo),PerfmIndxLib(stepNo),...
        hypParams,initVals,muGt_1,cGt_1);
    % update muGt and cGt before going to the next step
    muGt_1 = muGt;
    cGt_1 = cGt;
end

% calculate predicted mean and variance using traditional GP
mu = NaN*basisVecs;
sig = NaN*basisVecs;
for ii = 1:length(basisVecs)
    [mu(ii),sig(ii)] = rgp.calcPredMeanAndPredVar(basisVecs(ii),desPtSeq,...
        PerfmIndxLib,hypParams);
end

%% post-process
upperBound = muGt + diag(cGt);
lowerBound = muGt - diag(cGt);

%% plots
figure
hold on
grid on
% sampled points
scatter(desPtSeq*(180/pi),PerfmIndxLib,'r','filled')
% predicted mean from RGP
plot(basisVecs*(180/pi),muGt(:,end))
% predicted mean from traditional GP
plot(basisVecs*(180/pi),mu,'--k','LineWidth',2)
% bouunds from RGP
plot(basisVecs*(180/pi),upperBound,'-r','LineWidth',2)
plot(basisVecs*(180/pi),lowerBound,'-r','LineWidth',2)
% labels, legends, and title
xlabel('$\theta$')
ylabel('$J_{inst}$')
legend('f(x)','RGP','GP')
title('Mean')

figure
hold on
grid on
% prediction variance from RGP
plot(basisVecs*(180/pi),diag(cGt),'LineWidth',2)
% prediction variance from tradtional GP
plot(basisVecs*(180/pi),sig,'--k','LineWidth',2)
% labels, legends, and title
xlabel('$\theta$')
ylabel('$\sigma ^2 (\theta)$')
legend('RGP','GP')
title('Variance')


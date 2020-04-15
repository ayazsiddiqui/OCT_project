clear
clc
format compact

%% initialization
% construct an instance of the RGP class
rgp = RGP(1);

% set grid size
size_grid = 201;

% make basis vector space
basisVec = linspace(0,20,size_grid)*(pi/180);
% populate initial space with dummy values for the latent function
dummyVals = rand(size_grid,1);

% know values
desPtSeq = [0 0 1 2 3 4 15 5 7 6 8 12 8 8 8 8 8]*(pi/180);
PerfmIndxLib = [0.3 0.28 0.22 0.35 0.18 0.4 0.33 0.42 0.24 0.26 0.19 0.25 0.24 0.29 0.3 0.31 0.19]';

desPtSeq = desPtSeq(1:end-6);
PerfmIndxLib = PerfmIndxLib(1:end-6);

% set values of hyper parameters
hyp.funcVar = 1;         % variance of latent function
hyp.lengthScale = 0.01;  % length scale
hyp.noiseVar = 0.12 ;    % noise variance
hypParams = [hyp.funcVar;hyp.noiseVar;hyp.lengthScale];

% test function that calculates covariance and mean
[K,M] = rgp.buildCovMatAndMeanVec(basisVec,hypParams);

% test RGP
[initCovMat,initMeanVec] = ...
    rgp.buildCovMatAndMeanVec(basisVec,hypParams);
initInvCovMat = inv(initCovMat + (1e-5)*eye(size_grid));

cGt_11 = initCovMat;
muGt_11 = initMeanVec;

muGt_12 = muGt_11;
cGt_12 = cGt_11;

for stepNo = 1:length(desPtSeq)
    [muGt1,cGt1] = rgp.rgpRegression(desPtSeq(stepNo),PerfmIndxLib(stepNo),...
        hypParams,basisVec,initMeanVec,initInvCovMat,muGt_11,cGt_11);
    % update muGt and cGt before going to the next step
    muGt_11 = muGt1;
    cGt_11 = cGt1;
    
    %
%     [muGt2,cGt2,newbasisVec] = rgp.timeDependentRgpRegression(desPtSeq(stepNo),PerfmIndxLib(stepNo),...
%         hypParams,basisVec,muGt_12,cGt_12);
    
%     muGt_12 = muGt2;
%     cGt_12 = cGt2;

end

% calculate predicted mean and variance using traditional GP
mu = NaN*basisVec;
sig = NaN*basisVec;
for ii = 1:length(basisVec)
    [mu(ii),sig(ii)] = rgp.calcPredMeanAndPredVar(basisVec(ii),desPtSeq,...
        PerfmIndxLib,hypParams);
end

%% post-process
upperBound = muGt1 + diag(cGt1);
lowerBound = muGt1 - diag(cGt1);

%% plots
figure
hold on
grid on
% sampled points
scatter(desPtSeq*(180/pi),PerfmIndxLib,'r','filled')
% predicted mean from RGP
plot(basisVec*(180/pi),muGt1(:,end))
% predicted mean from traditional GP
plot(basisVec*(180/pi),mu,'--k','LineWidth',2)
% bouunds from RGP
plot(basisVec*(180/pi),upperBound,'-r','LineWidth',2)
plot(basisVec*(180/pi),lowerBound,'-r','LineWidth',2)
% labels, legends, and title
xlabel('$\theta$')
ylabel('$J_{inst}$')
legend('f(x)','RGP','GP')
title('Mean')

figure
hold on
grid on
% prediction variance from RGP
plot(basisVec*(180/pi),diag(cGt1),'LineWidth',2)
% prediction variance from tradtional GP
plot(basisVec*(180/pi),sig,'--k','LineWidth',2)
% labels, legends, and title
xlabel('$\theta$')
ylabel('$\sigma ^2 (\theta)$')
legend('RGP','GP')
title('Variance')


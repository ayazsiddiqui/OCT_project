close all
clear all
format compact
clc

% repeated random numbers
rng('default');
rng(1);

% original grd
nGridPt = 50;
x1 = linspace(0,5,nGridPt);
X1 = x1;
objF = @(x) 0.25*x.*sin(2*x).*exp(-0.5*x) - 0.125*x.^2.*exp(-x) + 0.12;

Z = objF(X1);
Zcol = Z(:);
nGrid = numel(Zcol);

XGrid = X1(:)';

% sampled data
nSampGrid = 10;
Fnoise = 0*0.025;
x1Sampled = linspace(0,5,nSampGrid);

X1Samp = x1Sampled;
ZSamp = objF(X1Samp);
ZSamp = ZSamp + 1*normrnd(0,Fnoise,size(ZSamp));

XSamp = X1Samp(:)';
ZSampCol = ZSamp(:);
nSamp = numel(ZSampCol);

% optimize hyper parameters
ini_sigma0 = 0.05;

% initial theta
a = 0;
b = 4;
ini_theta = ((b-a).*rand(size(XSamp,1),1) + a);
ini_hyperParam = [ini_sigma0;ini_theta];
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1+numel(ini_theta),1);
ub = Inf(1+numel(ini_theta),1);

sol_hyper = fmincon(@(hyperParam) logLogLikelihood(ZSampCol,XSamp,hyperParam),ini_hyperParam,A,b,Aeq,beq,lb,ub);

opt_sigma0 = sol_hyper(1);
opt_theta = sol_hyper(2:end);
% opt_sigma0 = 0.09;
% opt_theta = 1;
% Calculation of true predictive mean and variance from traditional GP
% modelling
trueCovMatrix = zeros(nSamp,nSamp);
for ii = 1:nSamp
    for jj = ii:nSamp
        trueCovMatrix(ii,jj) = covFuncEval(XSamp(:,ii),XSamp(:,jj),opt_sigma0,opt_theta);
    end
end
trueCovMatrix = trueCovMatrix + trueCovMatrix' - eye(size(trueCovMatrix)).*diag(trueCovMatrix);

% calculate
covRowVectorTrue = zeros(nSamp,nGrid);
for ii = 1:nGrid
    for jj = 1:nSamp
        covRowVectorTrue(jj,ii) = covFuncEval(XGrid(ii),XSamp(jj),opt_sigma0,opt_theta);
    end
end

noiseMat = 0.00025*eye(nSamp);
truePredMeanFunc = NaN(nGrid,1);
truePredVarFunc = NaN(nGrid,1);

for kk = 1:nGrid
    truePredMeanFunc(kk,1) = (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*ZSampCol(1:nSamp);
    truePredVarFunc(kk,1) = opt_sigma0 - (covRowVectorTrue(:,kk)'/(trueCovMatrix+ noiseMat))*covRowVectorTrue(:,kk);
end

upperLimitMean = truePredMeanFunc + 2*truePredVarFunc;
lowerLimitMean = truePredMeanFunc - 2*truePredVarFunc;

[minVal,minIndex] = min(truePredMeanFunc);

figure
% set(gca,'Xticklabel',[],'Yticklabel',[])
hold on
grid on
plot(X1,Z);
plot(X1(:),truePredMeanFunc,'r+')
plot(X1,upperLimitMean,'-.m')
plot(X1,lowerLimitMean,'-.m')

ps = scatter(X1Samp(:),ZSampCol,80,'filled','r');
legend(ps,{'Sampled data'});
hold off




%--------------------------------------------------------------------------




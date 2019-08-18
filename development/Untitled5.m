close all
clear all
format compact
clc

% original grd
nGridPt = 50;
x1 = linspace(-5,5,nGridPt);
x2 = linspace(-5,5,nGridPt);
X = [x1 x2];
objF = @(x1,x2) -(x1.^2 + x2.^2)./50 + 1;

[X1,X2] = meshgrid(x1,x2);
Z = objF(X1,X2);
Zcol = Z(:);
nGrid = numel(Zcol);

XGrid = [X1(:) X2(:)]';

% sampled data
nSampGrid = 10;
Fnoise = 0*0.025;
x1Sampled = linspace(-5,5,nSampGrid);
x2Sampled = linspace(-5,5,nSampGrid);

[X1Samp,X2Samp] = meshgrid(x1Sampled,x2Sampled);
ZSamp = objF(X1Samp,X2Samp);
ZSamp = ZSamp + 0*normrnd(0,Fnoise,size(ZSamp));

XSamp = [X1Samp(:) X2Samp(:)]';
ZSampCol = ZSamp(:);
nSamp = numel(ZSampCol);

% optimize hyper parameters
ini_sigma0 = 0.1;
a = 0.1;
b = 2;
ini_theta = (b-a).*rand(size(XSamp,1),1) + a;
ini_hyperParam = [ini_sigma0;ini_theta];

sol_hyper = fminunc(@(hyperParam) logLogLikelihood(ZSampCol,XSamp,hyperParam),ini_hyperParam);

opt_sigma0 = sol_hyper(1);
opt_theta = sol_hyper(2);

% Calculation of true predictive mean and variance from traditional GP
% modelling
trueCovMatrix = zeros(nSamp,nSamp);
for ii = 1:nSamp
    for jj = ii:nSamp
        trueCovMatrix(ii,jj) = covFuncEval(XSamp(:,ii),XSamp(:,jj),opt_theta,opt_sigma0);
    end
end
trueCovMatrix = trueCovMatrix + trueCovMatrix' - eye(size(trueCovMatrix)).*diag(trueCovMatrix);

% calculate
covRowVectorTrue = zeros(nSamp,nGrid);
for ii = 1:nGrid
    for jj = 1:nSamp
        covRowVectorTrue(jj,ii) = covFuncEval(XGrid(ii),XSamp(jj),opt_theta,opt_sigma0);
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
surf(X1,X2,Z);
view(-40,35);
% scatter3(X1(:),X2(:),truePredMeanFunc,'r+')
% surf(X1,X2,reshape(truePredMeanFunc,nGridPt,nGridPt));
% surf(X1,X2,reshape(upperLimitMean,nGridPt,nGridPt));
% surf(X1,X2,reshape(lowerLimitMean,nGridPt,nGridPt));
colorbar

ps = scatter3(X1Samp(:), X2Samp(:),ZSampCol,80,'filled','r');
legend(ps,{'Sampled data'});
hold off

figure
plot(1:nGrid,Zcol);
hold on
plot(1:nGrid,truePredMeanFunc);



%--------------------------------------------------------------------------




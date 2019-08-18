close all
clear all
format compact

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
Fnoise = 0.025;
x1Sampled = linspace(-5,5,nSampGrid);
x2Sampled = linspace(-5,5,nSampGrid);

[X1Samp,X2Samp] = meshgrid(x1Sampled,x2Sampled);
ZSamp = objF(X1Samp,X2Samp);
ZSamp = ZSamp + 0*normrnd(0,Fnoise,size(ZSamp));

XSamp = [X1Samp(:) X2Samp(:)]';
ZSampCol = ZSamp(:);
nSamp = numel(ZSampCol);

% hyper parameters
sigmaSq = 1;
hypPara = 0.5;

% ----------------------------------------------------------------------------
% Calculation of true predictive mean and variance from traditional GP
% modelling
trueCovMatrix = zeros(nSamp,nSamp);
for ii = 1:nSamp
    for jj = ii:nSamp
        trueCovMatrix(ii,jj) = covFuncEval(XSamp(:,ii),XSamp(:,jj),hypPara,sigmaSq);
    end
end
trueCovMatrix = trueCovMatrix + trueCovMatrix' - eye(size(trueCovMatrix)).*diag(trueCovMatrix);

% calculate
covRowVectorTrue = zeros(nSamp,nGrid);
for ii = 1:nGrid
    for jj = 1:nSamp
        covRowVectorTrue(jj,ii) = covFuncEval(XGrid(ii),XSamp(jj),hypPara,sigmaSq);
    end
end

noiseMat = 0.00025*eye(nSamp);
truePredMeanFunc = NaN(nGrid,1);
truePredVarFunc = NaN(nGrid,1);

for kk = 1:nGrid
    truePredMeanFunc(kk,1) = (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*ZSampCol(1:nSamp);
    truePredVarFunc(kk,1) = sigmaSq - (covRowVectorTrue(:,kk)'/(trueCovMatrix+ noiseMat))*covRowVectorTrue(:,kk);
end

upperLimitMean = truePredMeanFunc + 2*truePredVarFunc;
lowerLimitMean = truePredMeanFunc - 2*truePredVarFunc;

[minVal,minIndex] = min(truePredMeanFunc);

figure
% set(gca,'Xticklabel',[],'Yticklabel',[])
hold on
grid on
% surf(X1,X2,Z);
view(-40,35);
surf(X1,X2,reshape(truePredMeanFunc,nGridPt,nGridPt));
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
function [covFuncValue] = covFuncEval(designPt1,designPt2,theta,noiseVar)
covFuncValue = noiseVar*exp(-((designPt1-designPt2)'*(designPt1-designPt2))/(2*theta^2));
end



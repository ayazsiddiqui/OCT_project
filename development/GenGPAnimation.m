close all
clear all

x = [0:0.01:5]';
funcX = 0.25*x.*sin(2*x).*exp(-0.5*x) - 0.125*x.^2.*exp(-x) + 0.12;
size_grid = length(x);

sigmaSq = 0.125;
hypPara = 0.18;
xSampled = [0.25 4.26 3.75 0.95]';
funcXSampled = 0.25*xSampled.*sin(2*xSampled).*exp(-0.5*xSampled) - 0.125*xSampled.^2.*exp(-xSampled) + 0.12;
offlineSimEndSample = length(xSampled);
%--------------------------------------------------------------------------
% Calculation of true predictive mean and variance from traditional GP
% modelling

for ii = 1:offlineSimEndSample
    for jj = 1:offlineSimEndSample
        trueCovMatrix(ii,jj) = covFuncEval(xSampled(ii),xSampled(jj),hypPara,sigmaSq); 
    end 
end 

for ii = 1:size_grid
    for jj = 1:offlineSimEndSample
        covRowVectorTrue(jj,ii) = covFuncEval(x(ii),xSampled(jj),hypPara,sigmaSq);
    end 
end 

noiseMat = 0.009*eye(offlineSimEndSample);
for kk = 1:size_grid
%     trueMeanFunc(kk,1) = covRowVectorTrue(:,kk)'*inv(trueCovMatrix)*PerfmIndxLib(1:offlineSimEndSample);
%     trueCovFunc(kk,1) = 1 - covRowVectorTrue(:,kk)'*inv(trueCovMatrix)*covRowVectorTrue(:,kk);
    truePredMeanFunc(kk,1) = (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*funcXSampled(1:offlineSimEndSample);
    truePredVarFunc(kk,1) = sigmaSq - (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*covRowVectorTrue(:,kk);
end  

upperLimitMean = truePredMeanFunc + 2*truePredVarFunc;
lowerLimitMean = truePredMeanFunc - 2*truePredVarFunc;

figure
set(gca,'Xticklabel',[],'Yticklabel',[])
hold on 
grid on
plot(x,funcX,'--k')
plot(x,truePredMeanFunc)
plot(x,upperLimitMean,'-.m')
plot(x,lowerLimitMean,'-.m')
scatter(xSampled,funcXSampled,160,'filled','r')
legend('True func','Estimated mean','Upper Conf. int.','Lower Conf. int.','Sampled data')
hold off

%--------------------------------------------------------------------------    
function [covFuncValue] = covFuncEval(designPt1,designPt2,theta,noiseVar)
    covFuncValue = noiseVar*exp(-1/(2*theta^2)*(abs(designPt1-designPt2))^2);
end 
    
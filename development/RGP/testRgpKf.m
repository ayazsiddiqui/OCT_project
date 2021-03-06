clear
clc
format compact

%% generate wind using colored noise
% environment
hMax = 1500;
hMin = 100;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;
noTP = numel(heights);
% time in minutes
tVec = 0:2:6*60;
noTimeSteps = numel(tVec);
% time in seconds
timeInSec = 60*tVec;
% std deviation for wind data generation
stdDev = 0.9;
% hyper parameters
timeScale = 10;
heightScale = 200;
% generate data
windSpeedOut = meanFlow*(1 + genWindv2(heights,heightScale,tVec,timeScale,stdDev));

% p = polyfitn([repmat(heights,noTimeSteps,1), repmat(tVec,noTimeSteps,1)],...
%     windSpeedOut(:,1),3);

%% set up the RGP
% construct an instance of the RGP class
rgp = RGP(2);
% make basis vector space
basisVec = [heights(:)'; zeros(1,numel(heights))];
% set values of hyper parameters
hyp.funcVar = 1;         % variance of latent function
hyp.lengthScale = 0.01;  % length scale
hyp.timeScale = 5;       % time scale
hyp.noiseVariance = 0.12 ;    % noise variance
hyperParams = [hyp.funcVar;hyp.noiseVariance;hyp.lengthScale;hyp.timeScale];
% test log likelihood calculation
logP = rgp.calcMarginalLikelihood(basisVec,windSpeedOut(:,1),hyperParams)

% optimize hyper parameters by maxmizing the marginal likelihood
options = optimoptions('fmincon');
[optHyperParams,minLogP] = fmincon(@(hyperParams)...
    -rgp.calcMarginalLikelihood(basisVec,windSpeedOut(:,1),hyperParams),...
    hyperParams,[],[],[],[],hyperParams*0,[],[],options);


optHyperParams(3) = 200;
optHyperParams(4) = 10;
optHyperParams(2) = 1;


%% perform the actual recurssion
% number of test points at each step
nTest = 5;

% repeat the basis vector nRepeat times
nRepeat = 3;
newbasisVec = repmat(basisVec,1,nRepeat);

% initialize values of the covariance matrix and mean function
[initVals.initCovMat,initVals.initMeanVec] = ...
    rgp.buildCovMatAndMeanVec(newbasisVec,optHyperParams);

cGt_1 = initVals.initCovMat;
muGt_1 = initVals.initMeanVec;

% preallocate matrices
predMean = NaN(noTP,noTimeSteps);
upperBound = NaN(noTP,noTimeSteps);
lowerBound = NaN(noTP,noTimeSteps);
predVar =  NaN(noTP,noTimeSteps);

for stepNo = 1:noTimeSteps
    
    % randomly pick nTest points from the altitude basis vector
    tpIdx = randi([1 noTP],nTest,1);
    
    for jj = 1:nTest
        % update muGt and cGt before going to the next step
        [muGt,cGt,newbasisVec] = rgp.timeDependentRgpRegression([heights(tpIdx(jj));...
            tVec(stepNo)],windSpeedOut(tpIdx(jj),stepNo),...
            optHyperParams,newbasisVec,muGt_1,cGt_1);
        
    end
    
    muGt_1 = muGt;
    cGt_1 = cGt;
    
    predMean(:,stepNo) = mean(reshape(muGt,noTP,[]),2);
    predVar(:,stepNo) = mean(reshape(diag(cGt),noTP,[]),2);
    upperBound(:,stepNo) = predMean(:,stepNo) + predVar(:,stepNo);
    lowerBound(:,stepNo) = predMean(:,stepNo) - predVar(:,stepNo);
    
end

%% tradional GP estimates
% calculate predictive mean and variance
% mu = NaN(noTP,1);
% sig = NaN(noTP,1);
%
% for ii = 1:noTP
%     [mu(ii),sig(ii)] = rgp.calcPredMeanAndPredVar(basisVec(:,ii),basisVec,...
%         windSpeedOut(:,1),optHyperParams);
% end
% upperBound = mu + sig;
% lowerBound = mu - sig;

%% plot data
lwd = 1;

% set find plot limits
lB = min([lowerBound(:);windSpeedOut(:)]);
uB = max([upperBound(:);windSpeedOut(:)]);
plotRes = 5;

figure
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);

for ii = 1:noTimeSteps
    
    if ii == 1
        hold on
        grid on
        xlabel('Wind speed (m/s)');
        ylabel('Altitude (m)');
        
        %         xlim([lB-mod(lB,plotRes),uB-mod(uB,plotRes)+plotRes])
        xlim(20*[0 1.2])
        ylim([hMin hMax]);
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','color','k','-or','color','r','-or','color','b');
        delete(h);
        
    end
    
    plot(windSpeedOut(:,ii),heights,'k','linewidth',lwd);
    plot(predMean(:,ii),heights,'r','linewidth',lwd);
    plot(lowerBound(:,ii),heights,'b--','linewidth',lwd);
    plot(upperBound(:,ii),heights,'b--','linewidth',lwd);
    
    legend('True func','Pred mean','Bounds')
    txt = sprintf('Time = %0.2f min',tVec(ii));
    title(txt);
    
    
    
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

%%
% % % % video setting
video = VideoWriter('vid_Test1','Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 1;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)
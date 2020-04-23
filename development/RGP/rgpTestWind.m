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
% windSpeedOut = meanFlow*(1 + genWindv2(heights,heightScale,tVec,timeScale,stdDev));
windSpeedOut = genWindv2(heights,heightScale,tVec,timeScale,stdDev);
heights2 = repmat(heights,1,noTimeSteps);
tVec2 = repmat(tVec(:)',noTP,1);
dsgnPts = [heights2(:)'; tVec2(:)'];
dsgnFvals = windSpeedOut(:);


% translate and scale wind data
Flows = NaN(noTP,1,noTimeSteps);
for ii = 1:noTimeSteps
    Flows(:,:,ii) = windSpeedOut(:,ii);
end

simTime = 1*60*60;
Heights = timeseries(repmat(heights,1,1,2),[0 simTime]);
Flows = timeseries(Flows,timeInSec);


%% set up the RGP
% construct an instance of the RGP class
rgp = RGP(1);
% make basis vector space
basisVec = heights(:)';
% set values of hyper parameters
hyp.funcVar = 1;         % variance of latent function
hyp.lengthScale = 0.01;  % length scale
hyp.timeScale = 5;       % time scale
hyp.noiseVariance = 0.12 ;    % noise variance
hyperParams = [hyp.funcVar;hyp.noiseVariance;hyp.lengthScale];
% test log likelihood calculation
logP = rgp.calcMarginalLikelihood(basisVec,windSpeedOut(:,1),hyperParams);

% optimize hyper parameters by maxmizing the marginal likelihood
options = optimoptions('fmincon');
[optHyperParams,minLogP] = fmincon(@(hyperParams)...
    -rgp.calcMarginalLikelihood(basisVec,windSpeedOut(:,1),hyperParams),...
    hyperParams,[],[],[],[],hyperParams*0,[],[],options);


optHyperParams(2) = 0.2;
optHyperParams(3) = 200;
optHyperParams(4) = 10;

%% perform the actual recurssion
% number of test points at each step
nTest = 1;

% initialize values of the covariance matrix and mean function
[initCovMat,initMeanVec] = ...
    rgp.buildCovMatAndMeanVec(basisVec,optHyperParams);

initInvCovMat = inv(initCovMat + optHyperParams(2)*eye(noTP));

cGt_1 = initCovMat;
muGt_1 = initMeanVec;
hyperParams = optHyperParams;


%% plot data
keyboard

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
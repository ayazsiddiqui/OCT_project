clear
clc
format compact
close all

%% generate wind using colored noise
rng(1);

% environment
hMax = 1500;
hMin = 100;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;
noTP = numel(heights);
% time in minutes
timeStep = 0.05;
tVec = 0:timeStep:1*60;
noTimeSteps = numel(tVec);
% time in seconds
timeInSec = 60*tVec;
% std deviation for wind data generation
stdDev = 0.5;
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

%% set up the KFGP classdef
% construct an instance of the RGP class
gpkf = GPKF(1);
% set values of hyper parameters
hyperParams = [1 heightScale timeScale 0*0.0001]';
% % % optimize hyper parameters by maxmizing the marginal likelihood
% options = optimoptions('fmincon');
% [optHyperParams,minLogP] = fmincon(@(hyperParams)...
%     -gpkf.calcMarginalLikelihood(dsgnPts,dsgnFvals,hyperParams),...
%     hyperParams,[],[],[],[],[0;heightScale;timeScale;0.0],...
%     [Inf;heightScale;timeScale;Inf],[],options);

misCalc = (0.5-0).*rand(4,1) + 0;
optHyperParams = 1.*[1 heightScale timeScale 0*0.0001]';
hyperError = optHyperParams./hyperParams;
% optHyperParams = [38  heightScale timeScale 0.004]';
% % % test log likelihood calculation
% logP = gpkf.calcMarginalLikelihood(dsgnPts,dsgnFvals,optHyperParams);


%% do the actual gaussian process kalman filtering
% % % make domain vector
xDomain = heights(:)';
% % % set the measurable domain equal to the entire domain
xMeasure = xDomain;
% % % make a finer domain over which predictions are made
xPredict = linspace(heights(1),heights(end),1*numel(heights));
% % % order of approixation for SE kernel
Nn = 2;
% form the initialization matrices
initCons = gpkf.seGpkfInitialize(xDomain,optHyperParams(end-1),timeStep,Nn);
% % % set number of points visited per step
nVisit = 1;
% % % initialize parameters
noiseVar = optHyperParams(end);
Ks = gpkf.buildSpatialCovMat(xMeasure,optHyperParams(1),optHyperParams(2));

Ks_12 = chol(Ks,'upper');
Ks_12 = Ks_12 + triu(Ks_12,1)';
% Ks_12 = sqrtm(Ks);

ck_k = initCons.sig0Mat;
sk_k = initCons.s0;

% % % number of iterations
noIter = noTimeSteps;
% % % preallocate matrices
predMean = NaN(size(xPredict,2),noIter);
postVar =  NaN(size(xPredict,2),noIter);
stdDev =  NaN(size(xPredict,2),noIter);
upperBound = NaN(size(xPredict,2),noIter);
lowerBound = NaN(size(xPredict,2),noIter);
pointsVisited = NaN(nVisit,noIter);
fValAtPt = NaN(nVisit,noIter);



for ii = 1:noIter
    % % % visit said points
    if ii == 1
        visitIdx = sort(randperm(size(xMeasure,2),nVisit));
    else
        [~,sortIdx] =  max(postVar(:,ii-1));
        [~,visitIdx] = min((xMeasure - xPredict(sortIdx)).^2);
    end
    % % % extract visited values from xMeasure
    Mk = xMeasure(:,visitIdx);
    % % % extract wind speed at visited values
    yk = windSpeedOut(visitIdx,ii);
    % % % stepwise update of kalman state estimate and covariance
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
        gpkf.gpkfKalmanEstimation(xMeasure,sk_k,ck_k,Mk,yk,...
        Ks_12,initCons.Amat,initCons.Qmat,initCons.Hmat,noiseVar);
    % % % regression over a finer domain
    [predMean(:,ii),postVar(:,ii)] = gpkf.gpkfRegression(xDomain,xPredict,...
        F_t,sigF_t,Ks,optHyperParams);
    % % % remove real or imaginary parts lower than eps
    stdDev(:,ii) = sqrt(gpkf.removeEPS(postVar(:,ii),5));
    % % % upper bounds = mean + x*(standard deviation)
    upperBound(:,ii) = predMean(:,ii) + 1*stdDev(:,ii);
    % % %lower bounds = mean + x*(standard deviation)
    lowerBound(:,ii) = predMean(:,ii) - 1*stdDev(:,ii);
    % % % store points visited at the respective function value
    pointsVisited(:,ii) = Mk(:);
    fValAtPt(:,ii) = yk(:);
    % % % update previous step information
    sk_k = skp1_kp1;
    ck_k = ckp1_kp1;
    
end



%% plot data
% keyboard

lwd = 1;

% % % set find plot limits
lB = min([lowerBound(:);windSpeedOut(:)]);
uB = max([upperBound(:);windSpeedOut(:)]);
plotRes = 1;

figure(1)
set(gcf,'position',[1351 551 560 420]);
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);

for ii = 1:noTimeSteps
    
    if ii == 1
        hold on
        grid on
        xlabel('Wind speed (m/s)');
        ylabel('Altitude (m)');
        
        %                 xlim([lB-mod(lB,plotRes),uB-mod(uB,plotRes)+plotRes])
        xlim([-1 1])
        ylim([hMin hMax]);
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','color','k','-or','color','r',...
            '-or','color','b','-or','color','m');
        delete(h);
        
    end
    
    plot(windSpeedOut(:,ii),heights,'k','linewidth',lwd);
    plot(predMean(:,ii),xPredict,'r','linewidth',lwd);
    plot(lowerBound(:,ii),xPredict,'b--','linewidth',lwd);
    plot(upperBound(:,ii),xPredict,'b--','linewidth',lwd);
    plot(fValAtPt(:,ii),pointsVisited(:,ii),'mo','linewidth',lwd);
    
    legend('True func','Pred mean','Bounds')
    txt = sprintf('Time = %0.2f min',tVec(ii));
    title(txt);
    
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

%%
% % % % video setting
video = VideoWriter('vid_Test2','Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 5;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)



clear
clc
format compact
close all

%% random seed
rng(1);

%% make the space
% time step
timeStep = 0.2;
% final time
tf = 3;
% time vector
tVec = 0:timeStep:tf;
% spatial dicretization
xDomain = 1:1:100;

%% create true function sample space
trueFunc = NaN(numel(xDomain),numel(tVec));
for ii = 1:numel(tVec)
    trueFunc(:,ii) = rand(numel(xDomain),1);
end
    
%% given values
% time scale
timeScale = 1;
% length scale
lengthScale = 1;
% covariance amplitude
covAmp = 0.2;
% noise variance
noiseVar = 1;
% points visited at each iteration
nVisit = 80;

%% construct an instance of the GPKF class
gpkf = GPKF(1);
% initialize the GPKF
Nn = 6;
initCons = gpkf.squaredExponentialGpkfInitialize(xDomain,...
    timeScale,timeStep,Nn);
% build the spatial covariance matrix
Ks = gpkf.buildSpatialCovMat(xDomain,covAmp,lengthScale);
% principal square root
Ks_12 = chol(Ks,'upper');
Ks_12 = Ks_12 + triu(Ks_12,1)';
% Ks_12 = sqrtm(Ks);
% initial values
ck_k = initCons.sig0Mat;
sk_k = initCons.s0;


%% do the regression
GPKFpredMean = NaN(numel(xDomain),numel(tVec));
GPKFpostVar = NaN(numel(xDomain),numel(tVec));
GPpredMean = NaN(numel(xDomain),numel(tVec));
GPpostVar = NaN(numel(xDomain),numel(tVec));
pointsVisited = NaN(nVisit,numel(tVec));
fValAtPt = NaN(nVisit,numel(tVec));
percFitWrtGP  = NaN(1,numel(tVec));
GPFit  = NaN(1,numel(tVec));
GPKFfit  = NaN(1,numel(tVec));

for ii = 1:numel(tVec)
    visitIdx = randperm(numel(xDomain),nVisit);
    % % % extract visited values from xDomain
    Mk = xDomain(visitIdx);
    % % % extract wind speed at visited values
    yk = trueFunc(visitIdx,ii);
    % % % store points visited at the respective function value
    pointsVisited(:,ii) = Mk(:)';
    fValAtPt(:,ii) = yk(:);
    % % % kalman state estimates
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
        gpkf.gpkfKalmanEstimation(xDomain,sk_k,ck_k,Mk,yk,...
        Ks_12,initCons.Amat,initCons.Qmat,initCons.Hmat,noiseVar);
    % % % regression over a finer domain
    [GPKFpredMean(:,ii),GPKFpostVar(:,ii)] = gpkf.gpkfRegression(xDomain,...
        xDomain,F_t,sigF_t,Ks,[covAmp;lengthScale;timeScale;noiseVar]);
    % % % update previous step information
    sk_k = skp1_kp1;
    ck_k = ckp1_kp1;
    % % % traditional GP regression
    % data base of all locations visited
    tVisited = tVec(1:ii).*ones(nVisit,ii);
    xVisited = pointsVisited(:)';
    xVisited = [xVisited(1:nVisit*ii); tVisited(:)'];
    % function value at the above locations
    yVisited = fValAtPt(:);
    yVisited = yVisited(1:nVisit*ii);
    % the actual regression
    [GPpredMean(:,ii),GPpostVar(:,ii)] = gpkf.traditionalGpRegression(...
        xVisited,yVisited,xDomain,tVec(ii),...
        [covAmp;lengthScale;timeScale;noiseVar]);
    
    % % % calculate percentage fits
    % GPKF fit when compared to GP fit
    percFitWrtGP(ii) = (1 - ((norm(GPKFpredMean(:,ii)-GPpredMean(:,ii)))/...
        norm(GPpredMean(:,ii))));
    % GP fit when compared to true function
    GPFit(ii) = (1 - (norm(trueFunc(:,ii) - GPpredMean(:,ii)))/...
        norm(trueFunc(:,ii)));
    % GPKF fit when compared to true function
    GPKFfit(ii) = (1 - ((norm(trueFunc(:,ii) - GPKFpredMean(:,ii)))/...
        norm(trueFunc(:,ii))));
    disp(num2str(ii));
    
end


%% plots
% % % linewidths
lwd = 1;

figure(1)
set(gcf,'position',[1268 0 1.0*[560 2*420]]);

subplot(2,1,1)
plot(tVec,percFitWrtGP);
grid on
hold on
xlabel('Time (s)')
ylabel('Fit (\%)')

subplot(2,1,2)
plot(tVec,GPFit);
grid on
hold on
plot(tVec,GPKFfit);
xlabel('Time (s)')
ylabel('Fit (\%)')
legend('GP fit','GPKF fit')



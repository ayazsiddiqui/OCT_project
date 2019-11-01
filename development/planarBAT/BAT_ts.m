clear
clc
format compact

rng(1) % For reproducibility

%% parameters
% environment
rhoFluid = 1;
gravAcc = 9.81;

% body
mass = 2000;
BF = 1.5;
turbDia = 10;
Sref = 0.25*pi*turbDia^2;
dragCoeff = 0.5;

% tether
numNode = 4;
tetDia = 0.005;
tetRho = 1300;
tetCD = 0.4;
tetYoungs = 4e9;
tetZeta = 0.05;

%% initial conditions
NumSys = 1;

xOffset = 100;
operZ = 500;
iniPos = [xOffset*(0:NumSys-1); operZ*ones(1,NumSys)].*ones(2,NumSys);
gndNodePos = iniPos.*[ones(1,NumSys);zeros(1,NumSys)];
iniVel = repmat([0;0],[1 NumSys]);
iniTetherLength = 0.99*sqrt(sum((iniPos-gndNodePos).^2));

% controller
kp = 0.2;
kd = 4*kp;
tau = 20;

% winch
maxWinch = 0.5;

%% signals
%% environment
hMax = 1000;
hMin = 0;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;

% generate wind using colored noise
% time in minutes
tVec = 0:2:3*60;
% time in seconds
timeInSec = 60*tVec;
% std deviation for wind data generation
stdDev = 0.8;
% hyper parameters
timeScale = 5;
heightScale = 200;
% generate data
windSpeedOut = genWindv2(heights,heightScale,tVec,timeScale,stdDev);

nS = length(tVec);

% translate and scale wind data
for ii = 1:nS
    Flows(:,:,ii) = meanFlow*(1 + windSpeedOut(:,ii));
end

%
fn = 1;
figure(1);
for ii = 1:nS
    if ii == 1
        grid on
        hold on
        xlabel('Flow speed (m/s)')
        ylabel('Altitude (m)')
        xlim(max(ceil(abs(Flows)),[],'all')*[0 1]);
    else
        delete(pF)
    end
    
    pF = plot(Flows(:,:,ii),heights,'k');
    title(sprintf('Time = %0.1f min',tVec(ii)));
    %         waitforbuttonpress
end

%% train GP
gp = timeDepGaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.05;

% number of samples used to train GP
nTrain = 100;
trainDsgns = NaN(gp.noInputs,nTrain);
trainFval = NaN(nTrain,1);
iHeight = randi(numel(heights),1,nTrain);
iTime = randi(numel(timeInSec),1,nTrain);

for ii = 1:nTrain
    trainDsgns(:,ii) = [heights(iHeight(ii));timeInSec(iTime(ii))];
    trainFval(ii,1) = Flows(iHeight(ii),1,iTime(ii));
end
% rearrange values such that time is increasing: this isnt necessary
[trainDsgns,ia,~] = unique(trainDsgns(:,:).','rows');
trainDsgns = trainDsgns';
trainFval = trainFval(ia);

[~,I] = sort(trainDsgns(2,:));
trainDsgns = trainDsgns(:,I);
trainFval = trainFval(I,1);

%% simulate
simTime = 2*60*60;
heights = timeseries(repmat(heights,1,1,2),[0 simTime]);
Flows = timeseries(Flows,timeInSec);
% time interval between optimizations
optDt = 2*60;
% minimimum percentage improvement to increase optimization bounds
gamma = 0.01;
% optimization bounds increment factor
beta = 1.1;
% design limits
designLimits = [100 hMax];

% initial point
iniPt = [500;0];

% MPC parameters
maxIter = simTime/optDt;
predSteps = 3;
timeStep = optDt;
iniTauPerc = maxWinch*timeStep/range(designLimits);

%% simulate optimization
open_system('BAT_th','loadonly')
try
    set_param(gcs,'SimulationCommand','Update')
catch
end
sim('BAT_th')

parseLogsout
% resample
dt = 1;
tNew = 0:dt:simTime;
signals = fieldnames(tsc);

for ii = 1:length(signals)
    tscResampleOpt.(signals{ii}) = resample(tsc.(signals{ii}),tNew);
end

%% simulate baseline
altSP = 500;

open_system('baseline_th','loadonly')
try
    set_param(gcs,'SimulationCommand','Update')
catch
end
sim('baseline_th')

parseLogsout
% resample

for ii = 1:length(signals)
    tscResampleBase.(signals{ii}) = resample(tsc.(signals{ii}),tNew);
end
posDataOptBaseRn = tscResampleBase.allNodePos.Data;
posDataOptBaseRn = reshape(posDataOptBaseRn,2,numNode,[],length(tNew));

%%
posDataOptRn = tscResampleOpt.allNodePos.Data;
posDataOptRn = reshape(posDataOptRn,2,numNode,[],length(tNew));

plotMargin = 50;
[xmin,xmax] = bounds(squeeze(posDataOptRn(1,:,:,:)),'all');
[zmin,zmax] = bounds(squeeze(posDataOptRn(2,:,:,:)),'all');

bx = [xmin - mod(xmin,plotMargin) - plotMargin,...
    xmax - mod(xmax,plotMargin) + plotMargin];
bz = [zmin - mod(zmin,plotMargin) - plotMargin,...
    zmax - mod(zmax,plotMargin) + plotMargin];

lwd = 1;
boxWidth = 0.04;
boxHeight = 0.03;

fn = fn+1;
figure(fn)
set(gcf,'Position',[200 100 3*560 2*420]);

xAxLim = [-(plotMargin) max(abs(bx(:)))+(plotMargin)];
zAxLim = [hMin hMax];

nSubplots = [2,3];

for ii = 1:length(tNew)
    
    subplot(nSubplots(1),nSubplots(2),1);
    
    if ii == 1
        hold on
        grid on
        xlabel('X (m)');
        ylabel('Z (m)');
        xlim(xAxLim);
        ylim(zAxLim);
        axOpt = gca;
        axLocOpt = axOpt.Position;
        
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','color','k','-or','color','r','-or','color','b');
        delete(h);
    end
    
    for jj = 1:NumSys
        plot(posDataOptRn(1,:,jj,ii),posDataOptRn(2,:,jj,ii),'k-o','linewidth',lwd);
        annotation('textbox',...
            [axLocOpt(1)-(boxWidth/2)+(axLocOpt(3)*(posDataOptRn(1,end,jj,ii)-xAxLim(1))/(xAxLim(2)-xAxLim(1)))...
            axLocOpt(2)-(boxHeight/2)+(axLocOpt(4)*posDataOptRn(2,end,jj,ii)/(zAxLim(2)-zAxLim(1)))...
            boxWidth boxHeight],...
            'EdgeColor','k',...
            'BackgroundColor','w',...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'String',{['BAT ',num2str(jj)]},...
            'LineWidth',0.8);
        
        optSP = plot(xAxLim,tscResampleOpt.optAlt.Data(:,:,ii)*[1 1],'r');
        
        % plot flow on different axis
        subplot(nSubplots(1),nSubplots(2),2);
        
        if ii == 1
            grid on
            hold on
            xlabel('Flow speed (m/s)')
            ylabel('Altitude (m)')
            xlim(ceil(max(tscResampleOpt.flowVels.Data(:,:,:),[],'all'))*[0 1]);
            ylim([hMin hMax]);
        else
            h = findall(gca,'type','line','color','k','-or','color','r','-or','color','b');
            delete(h);
        end
        
        plot(squeeze(tscResampleOpt.flowVels.Data(:,:,ii)),...
            squeeze(tscResampleOpt.flowAlts.Data(:,:,ii)),'Color','b');
        
        % plot incident flow
        subplot(nSubplots(1),nSubplots(2),3);
        if ii == 1
            grid on
            hold on
            xlabel('Time (s)')
            ylabel('Incident flow (m/s)')
            xlim(simTime*[0 1]);
            ylim(ceil(max(tscResampleOpt.flowVels.Data(:,:,:),[],'all'))*[0 1]);
        end
        
        plot(squeeze(tscResampleOpt.incidentFlow.Time(ii))*[1 1],...
            squeeze(tscResampleOpt.incidentFlow.Data(:,:,ii))*[1 1],'.','Color','k');
        
    end
    
    % base line plots
    subplot(nSubplots(1),nSubplots(2),4);
    
    if ii == 1
        hold on
        grid on
        xlabel('X (m)');
        ylabel('Z (m)');
        xlim(xAxLim);
        ylim(zAxLim);
        axBase = gca;
        axLocBase = axBase.Position;
        
    else
        delete(findall(gca,'type','annotation'));
        h = findall(gca,'type','line','color','k','-or','color','r','-or','color','b');
        delete(h);
    end
    
    for jj = 1:NumSys
        plot(posDataOptBaseRn(1,:,jj,ii),posDataOptBaseRn(2,:,jj,ii),'k-o','linewidth',lwd);
        annotation('textbox',...
            [axLocBase(1)-(boxWidth/2)+(axLocBase(3)*(posDataOptBaseRn(1,end,jj,ii)-xAxLim(1))/(xAxLim(2)-xAxLim(1)))...
            axLocBase(2)-(boxHeight/2)+(axLocBase(4)*posDataOptBaseRn(2,end,jj,ii)/(zAxLim(2)-zAxLim(1)))...
            boxWidth boxHeight],...
            'EdgeColor','k',...
            'BackgroundColor','w',...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'String',{['BAT ',num2str(jj)]},...
            'LineWidth',0.8);
        
        optSP = plot(xAxLim,tscResampleBase.optAlt.Data(:,:,ii)*[1 1],'r');
        
        % plot flow on different axis
        subplot(nSubplots(1),nSubplots(2),5);
        
        if ii == 1
            grid on
            hold on
            xlabel('Flow speed (m/s)')
            ylabel('Altitude (m)')
            xlim(ceil(max(tscResampleBase.flowVels.Data(:,:,:),[],'all'))*[0 1]);
            ylim([hMin hMax]);
        else
            h = findall(gca,'type','line','color','k','-or','color','r','-or','color','b');
            delete(h);
        end
        
        plot(squeeze(tscResampleBase.flowVels.Data(:,:,ii)),...
            squeeze(tscResampleBase.flowAlts.Data(:,:,ii)),'Color','b');
        
        % plot incident flow
        subplot(nSubplots(1),nSubplots(2),6);
        if ii == 1
            grid on
            hold on
            xlabel('Time (s)')
            ylabel('Incident flow (m/s)')
            xlim(simTime*[0 1]);
            ylim(ceil(max(tscResampleBase.flowVels.Data(:,:,:),[],'all'))*[0 1]);
        end
        
        plot(squeeze(tscResampleBase.incidentFlow.Time(ii))*[1 1],...
            squeeze(tscResampleBase.incidentFlow.Data(:,:,ii))*[1 1],'.','Color','k');
        
    end
    
    % Title
    txt = sprintf('Time = %0.2f s',tNew(ii));
    supertitle(txt);
    
    F(ii) = getframe(gcf);
    
end

%%
% % % video setting
% video = VideoWriter('vid_Test1', 'Motion JPEG AVI');
% video.FrameRate = 30*1/dt;
% set(gca,'nextplot','replacechildren');
% 
% open(video)
% for i = 1:length(F)
%     writeVideo(video, F(i));
% end
% close(video)






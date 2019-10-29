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
hMin = 100;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;

% generate wind using colored noise
% time in minutes
tVec = 0:5:60;
% time in seconds
timeInSec = 60*tVec;
% std deviation for wind data generation
stdDev = 0.8;
% hyper parameters
timeScale = 30;
heightScale = 200;
% generate data
windSpeedOut = genWindv2(heights,heightScale,tVec,timeScale,stdDev);

nS = length(tVec);

% translate and scale wind data
for ii = 1:nS
    Flows(:,:,ii) = meanFlow*(1 + windSpeedOut(:,ii));
end

% 
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
%     waitforbuttonpress
end

%% train GP
gp = timeDepGaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.05;

tVals = reshape(transpose(timeInSec'.*ones(length(timeInSec),length(heights))),[],1)';
trainDsgns = [repmat(heights',1,nS);tVals];
trainFval = Flows(:);

%% simulate
simTime = 20*60;
heights = timeseries(repmat(heights,1,1,2),[0 simTime]);
Flows = timeseries(Flows,timeInSec);
% time interval between optimizations
optDt = 1.5*60;
% minimimum percentage improvement to increase optimization bounds
gamma = 0.01;
% optimization bounds increment factor
beta = 1.1;
% design limits
designLimits = [hMin hMax];

% initial point
iniPt = [500;tVals(end)];
iniTauPerc = 0.1;

% MPC parameters
maxIter = simTime/optDt;
predSteps = 5;
timeStep = optDt;

% simulate
open_system('BAT_th','loadonly')
try
    set_param(gcs,'SimulationCommand','Update')
catch
end
sim('BAT_th')

%% postprocess
parseLogsout

% resample
dt = 1;
tNew = 0:dt:simTime;
signals = fieldnames(tsc);

for ii = 1:length(signals)
    tscResample.(signals{ii}) = resample(tsc.(signals{ii}),tNew);
end

%%
posData = tscResample.allNodePos.Data;
posData = reshape(posData,2,numNode,[],length(tNew));

plotMargin = 50;
[xmin,xmax] = bounds(squeeze(posData(1,:,:,:)),'all');
[zmin,zmax] = bounds(squeeze(posData(2,:,:,:)),'all');

bx = [xmin - mod(xmin,plotMargin) - plotMargin,...
    xmax - mod(xmax,plotMargin) + plotMargin];
bz = [zmin - mod(zmin,plotMargin) - plotMargin,...
    zmax - mod(zmax,plotMargin) + plotMargin];

lwd = 1;
boxWidth = 0.09;
boxHeight = 0.05;

figure(1)
xAxLim = [-(plotMargin) max(abs(bx(:)))+(plotMargin)];
zAxLim = [0 max(bz(:)) + (plotMargin)];

for ii = 1:length(tNew)
    figure(1)
    
    if ii == 1
        hold on
        grid on
        xlabel('X (m)');
        ylabel('Z (m)');
        xlim(xAxLim);
        ylim(zAxLim);
        ax1 = gca;
        axLOc = ax1.Position;
        
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','color','k','-or','color','r','-or','color','b');
        delete(h);
    end
    
    for jj = 1:NumSys
        plot(posData(1,:,jj,ii),posData(2,:,jj,ii),'k-o','linewidth',lwd);
        annotation('textbox',...
            [axLOc(1)-(boxWidth/2)+(axLOc(3)*(posData(1,end,jj,ii)-xAxLim(1))/(xAxLim(2)-xAxLim(1)))...
            axLOc(2)-(boxHeight/2)+(axLOc(4)*posData(2,end,jj,ii)/(zAxLim(2)-zAxLim(1)))...
            boxWidth boxHeight],...
            'EdgeColor','k',...
            'BackgroundColor','w',...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'String',{['BAT ',num2str(jj)]},...
            'LineWidth',0.8);
        
        optSP = plot(xAxLim,tscResample.optAlt.Data(:,:,ii)*[1 1],'r');
        
        % plot flow on different axis
        figure(2)
        plot(squeeze(tscResample.flowVels.Data(:,:,ii)),...
            squeeze(tscResample.flowAlts.Data(:,:,ii)),'Color','b')

    end
    
    title(['Time = ',sprintf('%0.2f', tNew(ii)),' s'])
    F(ii) = getframe(gcf);
    
end
hold off

%%
% % % video setting
% video = VideoWriter('vid_Test', 'Motion JPEG AVI');
% video.FrameRate = 30*1/dt;
% set(gca,'nextplot','replacechildren');
% 
% open(video)
% for i = 1:length(F)
%     writeVideo(video, F(i));
% end
% close(video)






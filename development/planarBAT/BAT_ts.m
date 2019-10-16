clear
clc
format compact

rng(1) % For reproducibility

%% parameters
% environment
rhoFluid = 1;
gravAcc = 9.81;

% body
mass = 5000;
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

% initial conditions
iniPos = [0;500];
iniVel = [0;0];
iniTetherLength = 0.99*norm(iniPos);

% controller
kp = 2;
kd = 4*kp;
tau = 0.5;

% winch
maxWinch = 0.5;

%% signals
Vflow = [5;0];
SP = 500;

%% environment
heights = [0:25:1000]';
meanFlow = 10;


Flows = normrnd(meanFlow,0.1*meanFlow,size(heights));

%% simulate
simTime = 300;
heights = timeseries(repmat(heights,1,1,2),[0 simTime]);
Flows = timeseries(repmat(Flows,1,1,2),[0 simTime]);


sim('BAT_th')

%% postprocess
parseLogsout

% resample
dt = 1/10;
tNew = 0:dt:simTime;
signals = fieldnames(tsc);

for ii = 1:length(signals)
    tscResample.(signals{ii}) = resample(tsc.(signals{ii}),tNew);
end

%%
posData = tscResample.allNodePos.Data;

plotMargin = 50;
[xmin,xmax] = bounds(squeeze(posData(1,:,:)),'all');
[zmin,zmax] = bounds(squeeze(posData(2,:,:)),'all');

bx = [xmin - mod(xmin,plotMargin) - plotMargin,...
    xmax - mod(xmax,plotMargin) + plotMargin];
bz = [zmin - mod(zmin,plotMargin) - plotMargin,...
    zmax - mod(zmax,plotMargin) + plotMargin];

lwd = 1;

for ii = 1:length(tNew)
    
    
    if ii == 1
        hold on
        grid on
        xlabel('X (m)');
        ylabel('Z (m)');
        xlim([-max(abs(bx(:)))-(plotMargin) max(abs(bx(:)))+(plotMargin)]);
        ylim([0 max(bz(:)) + (plotMargin)]);
        
    else
        delete(p1);
    end
    
    p1 = plot(posData(1,:,ii),posData(2,:,ii),'k-o','linewidth',lwd);
    title(['Time = ',sprintf('%0.2f', tNew(ii)),' s'])
    
    F(ii) = getframe(gcf);
    
end

%%
% % % video setting
% video = VideoWriter('vid_Test', 'Motion JPEG AVI');
% video.FrameRate = 1*1/dt;
% set(gca,'nextplot','replacechildren');
% 
% open(video)
% for i = 1:length(F)
%     writeVideo(video, F(i));
% end
% close(video)






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

%% initial conditions
NumSys = 1;

xOffset = 100;
operZ = 500;
iniPos = [xOffset*(0:NumSys-1); operZ*ones(1,NumSys)].*ones(2,NumSys);
gndNodePos = iniPos.*[ones(1,NumSys);zeros(1,NumSys)];
iniVel = repmat([0;0],[1 NumSys]);
iniTetherLength = 0.99*sqrt(sum((iniPos-gndNodePos).^2));

% controller
kp = 2;
kd = 4*kp;
tau = 0.5;

% winch
maxWinch = 0.5;

%% signals
%% environment
hMax = 1000;
hMin = 0;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;

% fit a polynomial to dummy data
ht = hMin:200:hMax;
ht = ht(:);
ft2 = [5;5.5;6.5;7;6;5.5];

pt = polyfit(ht,ft2,3);

Flows = polyval(pt,heights);

% simTime in minutes
simTime = 30;
FlowInt = 60*(0:5:simTime);
nS = length(FlowInt);

for ii = 1:nS
    
    for jj = 1:length(pt)
        rdNum = rand;
        if rdNum < 0.33
            mut = -1;
        elseif rdNum >= 0.33 && rdNum < 0.66
            mut = 0;
        else
            mut = 1;
        end
        pt(jj) = pt(jj)*(1 + mut*0.1);
    end
    Flows(:,:,ii) = polyval(pt,heights);
end

% train GP
gp = timeDepGaussianProcess(2,'kernel','squaredExponential','acquisitionFunction','upperConfidenceBound');
if strcmpi(class(gp.acquisitionFunction),'acquisitionFunctions.upperConfidenceBound')
    gp.acquisitionFunction.explorationFactor = 1;
end
gp.kernel.noiseVariance = 1*0.05;

tVals = reshape(transpose(FlowInt'.*ones(length(FlowInt),length(heights))),[],1)';
trainDsgns = [repmat(heights',1,nS);tVals];
trainFval = Flows(:);

gamma = 0.01;
beta = 1.1;
designLimits = [hMin hMax];

iniPt = [500;tVals(end)];
iniTauPerc = 0.1;
maxIter = 5;
predSteps = 5;
timeStep = 100;


%% simulate
simTime = 20*60;
heights = timeseries(repmat(heights,1,1,2),[0 simTime]);
Flows = timeseries(Flows,FlowInt);
optDt = 5*60;

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
    
    if ii == 1
        hold on
        grid on
        xlabel('X (m)');
        ylabel('Z (m)');
        xlim(xAxLim);
        ylim(zAxLim);
        AxesHandle=findobj(gcf,'Type','axes');
        pt1 = get(AxesHandle,{'Position'});
        axLOc = pt1{:};
        
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','color','k');
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

    end
    
    title(['Time = ',sprintf('%0.2f', tNew(ii)),' s'])
    F(ii) = getframe(gcf);
    
end
hold off

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






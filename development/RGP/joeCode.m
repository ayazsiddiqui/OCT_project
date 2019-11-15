clear
clc
format compact

%% objective functions
% % % Park example 1
objF = @(X)-((X(1,:).^2 + X(2,:).^2)./50) + 1;
% % % Park example 2
% objF = @(X) 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
%     +0.75*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
% % https://www.hindawi.com/journals/mpe/2013/948303/ example
% objF = @(X) exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
%     exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
%     2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
%     1.5.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));

rngSeed = 46;
rng(rngSeed);

% test designs
xMin = -5; xMax = 5;
designLimits = [xMin*[1;1],xMax*[1;1]];

iniPt = ((xMax-xMin).*rand(2,1) + xMin);

[optDsgn,maxF] = multistartBFGS(objF,iniPt,designLimits(:,1),designLimits(:,2)...
    ,'nStarts',25);

%% post process
% make contour
nSamp = 50;
x = linspace(xMin,xMax,nSamp);
[x1s,x2s] = meshgrid(x,x);
fVal = objF([x1s(:)';x2s(:)']);
fVal = reshape(fVal,nSamp,[]);

figure(1)
contourf(x1s,x2s,fVal)
colorbar
hold on
plot(optDsgn(1,:),optDsgn(2,:),'-rs',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',10)






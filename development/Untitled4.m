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

noiseSigma = 0.009;

% train gaussian process
gprMdl = fitrgp(xSampled,funcXSampled,'Basis','constant'...
    ,'KernelFunction','squaredexponential','Sigma',noiseSigma,...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));

Loss = loss(gprMdl,x,funcX);

ypred = predict(gprMdl,xSampled);

%
figure
set(gca,'Xticklabel',[],'Yticklabel',[])
hold on 
grid on
plot(x,funcX,'--k')
plot(xSampled,ypred,'*')
% plot(x,upperLimitMean,'-.m')
% plot(x,lowerLimitMean,'-.m')
% scatter(xSampled,funcXSampled,160,'filled','r')
% legend('True func','Estimated mean','Upper Conf. int.','Lower Conf. int.','Sampled data')
% hold off


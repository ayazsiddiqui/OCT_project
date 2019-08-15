clear
clc
format compact

% load data file
load(fullfile(matlabroot,'examples','stats','gprdata2.mat'));

% fit gaussain regression process
gprMdl1 = fitrgp(x,y,'KernelFunction','squaredexponential');

% Find hyperparameters that minimize five-fold cross-validation loss by using automatic hyperparameter optimization.
rng default
gprMdl2 = fitrgp(x,y,'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'));

% predict values
ypred1 = resubPredict(gprMdl1);
ypred2 = resubPredict(gprMdl2);

% plot
figure();
plot(x,y,'r.');
hold on
plot(x,ypred1,'b');
plot(x,ypred2,'k','LineWidth',2);
xlabel('x');
ylabel('y');
legend({'data','Initial Fit','Optimized Fit'},'Location','Best');
title('Impact of Optimization');
hold off
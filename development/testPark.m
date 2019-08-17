clear
clc
format compact

% example from park

x10 = linspace(-5,5,10);
x20 = linspace(-5,5,10);

[x1,x2] = meshgrid(x10,x20);

Xtrain = [x1(:), x2(:)];

[objTbl,ytrain] = objectiveFunction(Xtrain);

figure(1)
subplot(1,2,1);
contourf(x1,x2,reshape(ytrain,numel(x10),numel(x10)));
colorbar

options = {'FitMethod','None','BasisFunction','constant','KernelFunction','squaredexponential'};

% train gaussian process
gprMdl = fitrgp(Xtrain,ytrain,'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','squaredexponential');

% test the gaussian process
x10_test = linspace(-5,5,50);
x20_test = linspace(-5,5,50);

[x1_test,x2_test] = meshgrid(x10_test,x20_test);

Xtest = [x1_test(:), x2_test(:)];
[objTbl2,ytest] = objectiveFunction(Xtest);

Loss = loss(gprMdl,Xtest,ytest);

ypred = predict(gprMdl,Xtest);

%% plot prediceted contour
figure(1)
subplot(1,2,2);
contourf(x1_test,x2_test,reshape(ypred,numel(x10_test),numel(x20_test)));
colorbar

%% optimization phase

% % Call bayesopt, capturing x and y in the objective function
% bo = bayesopt(@(T)objFcn(T,Xtrain,ytrain), kernel)
% % Define an objective function
% function Loss = objFcn(Vars, x, y)
% m = fitrgp(x, y, 'KernelFunction', char(Vars.KernelFunction), 'KFold', 5);
% Loss = kfoldLoss(m);
% end

clear
clc
format compact

load(fullfile(matlabroot,'examples','stats','gprdata2.mat'))
 
% Define kernel variable
kernel = optimizableVariable('KernelFunction',{'exponential','squaredexponential','matern32','matern52',...
'rationalquadratic','ardexponential','ardsquaredexponential','ardmatern32','ardmatern52','ardrationalquadratic'},...
'Type','categorical')

% Call bayesopt, capturing x and y in the objective function
bo = bayesopt(@(T)objFcn(T,x,y), kernel)
% Define an objective function
function Loss = objFcn(Vars, x, y)
m = fitrgp(x, y, 'KernelFunction', char(Vars.KernelFunction), 'KFold', 5);
Loss = kfoldLoss(m);
end
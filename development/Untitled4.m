clear
clc
format compact

load(fullfile(matlabroot,'examples','stats','gprdata.mat'));

% vInitialize length scales of the kernel function at 10
% and signal and noise standard deviations at the standard deviation of the response
sigma0 = std(ytrain);
sigmaF0 = sigma0;
d = size(Xtrain,2);
sigmaM0 = 10*ones(d,1);

% Fit the GPR model using the initial kernel parameter values. 
% Standardize the predictors in the training data. Use the exact fitting and prediction methods

gprMdl = fitrgp(Xtrain,ytrain,'Basis','constant','FitMethod','exact',...
'PredictMethod','exact','KernelFunction','ardsquaredexponential',...
'KernelParameters',[sigmaM0;sigmaF0],'Sigma',sigma0,'Standardize',1);

% Plot the log of learned length scales
figure()
plot((1:d)',log(sigmaM),'ro-');
xlabel('Length scale number');
ylabel('Log of length scale');


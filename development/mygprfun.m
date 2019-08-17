function [objective] = mygprfun(x,y)
gprMdl = fitrgp(x,y,'Basis','constant','FitMethod','exact',...
    'PredictMethod','exact','KernelFunction','squaredexponential');

objective = kfoldLoss(crossval(gprMdl));

end


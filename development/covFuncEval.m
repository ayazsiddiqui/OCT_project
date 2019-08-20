function [covFuncValue] = covFuncEval(designPt1,designPt2,sigma0,sigmaE,theta)

diagLamb = eye(numel(theta))./(theta.^2);

covFuncValue = sigma0*exp(-0.5*((designPt1-designPt2)'*diagLamb*(designPt1-designPt2))) + sigmaE*eye(numel(theta));

end
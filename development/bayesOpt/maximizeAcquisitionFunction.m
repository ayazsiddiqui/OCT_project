function [val,aFmax] = maximizeAcquisitionFunction(bestFval,trainDsgn,trainCovMat,trainFval,initialPt,bounds,noiseVar,hyper)

A = []; b = [];
Aeq = []; beq = [];

lb = bounds(:,1)';
ub = bounds(:,2)';
nonlcon = [];
options  = optimoptions('fmincon','Display','off');
[val,aFmax] = fmincon(@(optDsgn) ...
    -calcAcquisitionFunctionVal(optDsgn,bestFval,trainDsgn,trainCovMat,trainFval,noiseVar,hyper),...
    initialPt,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

function val = calcAcquisitionFunctionVal(postDsgn,bestFval,trainDsgn,trainCovMat,trainFval,noiseVar,hyper)

[predMean,predVar] = calcPredictiveMeanAndVariance(postDsgn,trainDsgn,trainCovMat,trainFval,noiseVar,hyper);

% http://krasserm.github.io/2018/03/21/bayesian-optimization/
fBest = bestFval;

stdDev = sqrt(predVar);
pd = makedist('Normal','mu',0,'sigma',1);
gm = gmdistribution(0,1);

if stdDev>0
    Z = (predMean-fBest)/stdDev;
    val = 1*(((predMean-fBest)*cdf(pd,Z)) + stdDev*pdf(gm,Z));
else
    val = 0;
end
end

% calculate predictive mean and variance
function [predMean,predVar] = calcPredictiveMeanAndVariance(postDsgn,trainDsgn,trainCovMat,trainFval,noiseVar,hyper)

CovMatVecTranspose = buildCovarianceMatrix(trainDsgn,postDsgn,hyper(1),noiseVar,hyper(2:end));

% mean and variance
nPost = size(postDsgn,2);

muD = NaN(nPost,1);
Var = NaN(nPost,1);
for ii = 1:nPost
    muD(ii,1) = (CovMatVecTranspose(:,ii)'/trainCovMat)*trainFval;
    Var(ii,1) = buildCovarianceMatrix(postDsgn(:,ii),postDsgn(:,ii),hyper(1),noiseVar,hyper(2:end))...
        - (CovMatVecTranspose(:,ii)'/trainCovMat)*CovMatVecTranspose(:,ii);
end

predMean = muD;
predVar = Var;

end
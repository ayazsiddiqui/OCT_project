function EI = calcExpectedImprovement(testDsgns,testY,optHyper,dsgnPt)

fBest = max(testY);

[predictedMean,predictedVariance] = posteriorCalc(testDsgns,testY,optHyper,dsgnPt);

%http://krasserm.github.io/2018/03/21/bayesian-optimization/
stdDev = sqrt(predictedVariance);
pd = makedist('Normal','mu',predictedMean,'sigma',stdDev);
gm = gmdistribution(predictedMean,stdDev);

if stdDev > 0
    Z = (predictedMean - fBest)/stdDev;
    EI = -1*((predictedMean - fBest)*cdf(pd,Z) + stdDev*pdf(gm,Z));
else
    EI = 0;
end

end
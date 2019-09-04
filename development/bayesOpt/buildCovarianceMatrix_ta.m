clear
clc

dsgnSet1 = rand(2,8);
dsgnSet2 = rand(2,1);
dsgn1Fval = rand(size(dsgnSet1,2),1);


covarianceAmp = 0.05;
noiseVariance = 0.1;
lengthScale = [1;1];

covMat = buildCovarianceMatrix(dsgnSet1,dsgnSet2,covarianceAmp,noiseVariance,lengthScale)
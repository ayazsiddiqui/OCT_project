clear
clc

dsgnSet1 = rand(2,8);
dsgnSet2 = rand(2,1);

covarianceAmp = 0.05;
noiseVariance = 0.1;
lengthScale = [1;1];

covMat = buildCovarianceMatrix(dsgnSet1,dsgnSet2,covarianceAmp,noiseVariance,lengthScale)

u = linspace(-5,5,100).*[1;1];

for ii = 1:size(u,2)
fn(ii) = ((u(1,ii)^2 + u(2,ii)^2)/50) + 1;

end
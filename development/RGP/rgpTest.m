clear
clc
format compact

%% test RGP from Huber's journal paper
x1 = [4 5 3];
fVal = [5; 8; 2];


lScale = 2.5;
covAmp = 1;
noiseVar = 1e-3;

K = buildTraingCovMat(x1,[],[],covAmp,lScale,noiseVar);
invK = inv(K);

% test point
x2 = 4;

% predictive mean calc
K_s = buildCovVec(x1,x2,covAmp,lScale);
muD = calculatePredMean(K_s,invK,fVal);

% covariance with itself
K_ss = cov(x2,x2,covAmp,lScale) + noiseVar;
sigmaSq = calculatePredVariance(K_ss,K_s,invK)

%% functions

%%%% covariance function
function val = cov(x1,x2,covAmp,lengthScl)
val = covAmp*(exp(-0.5*((x1-x2)'*(eye(numel(lengthScl))./(lengthScl.^2))*(x1-x2))));
end

%%%% build covariance matrix
function val = buildTraingCovMat(oldDsgns,oldMat,newDsgn,covAmp,lengthScl,noiseVar)

numOldDes = size(oldDsgns,2);

if isempty(oldMat)
    oldMat = zeros(numOldDes,numOldDes);
    
    for ii = 1:numOldDes
        for jj = ii:numOldDes
            oldMat(ii,jj) = cov(oldDsgns(:,ii),oldDsgns(:,jj),covAmp,lengthScl);
        end
    end
    oldMat = oldMat + triu(oldMat,1)' + eye(numOldDes)*noiseVar;
end

if isempty(newDsgn)
    val = oldMat;
else
    newDsgnCov = zeros(1,numOldDes);
    for ii = 1:numOldDes
        newDsgnCov(ii) = cov(oldDsgns(:,ii),newDsgn(:),covAmp,lengthScl);
    end
    
    covNewDsgn = cov(newDsgn,newDsgn,covAmp,lengthScl) + noiseVar;
    val = [oldMat newDsgnCov'; newDsgnCov covNewDsgn];
    
end

end

%%%% build covariance vector
function val = buildCovVec(oldDsgns,newDsgn,covAmp,lengthScl)

vecLengh = size(oldDsgns,2);
val = zeros(vecLengh,1);

for ii = 1:vecLengh
    val(ii) = cov(newDsgn(:),oldDsgns(:,ii),covAmp,lengthScl);
end
end

%%%% calcuate predicive mean
function val = calculatePredMean(covVec,invKmat,fVal)
val = covVec'*invKmat*fVal(:);
end

%%%% calcuate predicive variance
function val = calculatePredVariance(covItself,covVec,invKmat)
val = covItself - covVec'*invKmat*covVec;
end



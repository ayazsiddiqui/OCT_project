function [muD,Var] = posteriorCalc(testDsgns,testY,hyperParam,postDsgns)

CovMatVecTranspose = buildCovMat(testDsgns,postDsgns,'covAmplitude',hyperParam(1),'noiseVariance',hyperParam(2),'lengthScale',hyperParam(3:end));
testCovMat = buildCovMat(testDsgns,testDsgns,'covAmplitude',hyperParam(1),'noiseVariance',hyperParam(2),'lengthScale',hyperParam(3:end));

% CovMatVecTranspose
% muD = (CovMatVecTranspose'/testCovMat)*testY';
% Var = buildCovMat(postDsgns,postDsgns,'covAmplitude',optHyper(1),'noiseVariance',optHyper(2),'lengthScale',optHyper(3:end));

%% mean and variance
nPost = size(postDsgns,2);

muD = NaN(nPost,1);
Var = NaN(nPost,1);
for ii = 1:nPost
    muD(ii,1) = (CovMatVecTranspose(:,ii)'/testCovMat)*testY;
    Var(ii,1) = buildCovMat(postDsgns(:,ii),postDsgns(:,ii),'covAmplitude',hyperParam(1),'noiseVariance',...
        hyperParam(2),'lengthScale',hyperParam(3:end)) - (CovMatVecTranspose(:,ii)'/testCovMat)*CovMatVecTranspose(:,ii);
end

end

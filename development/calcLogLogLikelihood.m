function op = calcLogLogLikelihood(y,X,hyperParam)

covAmp = hyperParam(1);
noiseVar = hyperParam(2);
lengthScl = hyperParam(3:end);

Kmat = buildCovMat(X,X,'covAmplitude',covAmp,'noiseVariance',noiseVar,'lengthScale',lengthScl);

op = -1*(-0.5*(y'/Kmat*y) - 0.5*log(det(Kmat)));
% - 0.5*numel(y)*log(2*pi));

end
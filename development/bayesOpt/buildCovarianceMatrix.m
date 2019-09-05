%% other functions
% build covariance mat
function val = buildCovarianceMatrix(dsgnSet1,dsgnSet2,covarianceAmp,noiseVariance,lengthScale)

% function that builds covariance matrix between design sets
% inputs datapoints, covariance amplitude, noise variance and length scale

% length scale quantifies the relevance of the components in the input
% vector for predicting the outputs

% covariance formula
cov = @(x1,x2,covAmp,noiseVar,lengthScl)...
    covAmp*(exp(-0.5*((x1-x2)'*(eye(numel(lengthScl))./(lengthScl.^2))*(x1-x2))));

% convariance matrix size
covMatSize = [size(dsgnSet1,2), size(dsgnSet2,2)];
covMat = zeros(covMatSize);

if isequal(dsgnSet1,dsgnSet2)
    for ii = 1:covMatSize(1)
        for jj = ii:covMatSize(2)
            covMat(ii,jj) = cov(dsgnSet1(:,ii),dsgnSet2(:,jj),...
                covarianceAmp,noiseVariance,lengthScale);
            
        end
    end
    covMat = covMat + covMat' - eye(length(covMat)).*diag(covMat - noiseVariance);
    
else
    for ii = 1:covMatSize(1)
        for jj = 1:covMatSize(2)
            covMat(ii,jj) = cov(dsgnSet1(:,ii),dsgnSet2(:,jj),...
                covarianceAmp,noiseVariance,lengthScale);
        end
    end
end

val = covMat;

end

function covMat = buildCovMat(dsgnSet1,dsgnSet2,varargin)

% function that builds covariance matrix between design sets
% inputs datapoints, covariance amplitude, noise variance and length scale

% length scale quantifies the relevance of the components in the input
% vector for predicting the outputs

% parse input
p = inputParser;
addRequired(p,'dsgnSet1',@isnumeric);
addRequired(p,'dsgnSet2',@isnumeric);
addParameter(p,'covAmplitude',0,@isnumeric);
addParameter(p,'noiseVariance',0,@isnumeric);
addParameter(p,'lengthScale',1,@isnumeric);

parse(p,dsgnSet1,dsgnSet2,varargin{:});

if size(p.Results.lengthScale,1) ~= size(p.Results.dsgnSet1,1) || ...
        size(p.Results.lengthScale,1) ~= size(p.Results.dsgnSet2,1)
    error('Number of elements in length scale not equal to number of input variables')
end

% convariance matrix size
covMatSize = [size(dsgnSet1,2), size(dsgnSet2,2)];
covMat = zeros(covMatSize);

% choose kernel
kernel = 1;
switch kernel
    case 1
        % squared exponential
        cov = @(x1,x2,covAmp,noiseVar,lengthScl)...
            covAmp*(exp(-0.5*((x1-x2)'*(eye(numel(lengthScl))./(lengthScl.^2))*(x1-x2)))) + noiseVar*eye(numel(lengthScl));
end

for ii = 1:covMatSize(1)
    for jj = 1:covMatSize(2)
        covMat(ii,jj) = cov(dsgnSet1(:,ii),dsgnSet2(:,jj),...
            p.Results.covAmplitude,p.Results.noiseVariance,p.Results.lengthScale);
        
    end
end

end

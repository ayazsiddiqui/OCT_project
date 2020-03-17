classdef RGP
    %RGP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        noInputs
    end
    
    properties
        noBasisVec
    end
    
    
    %% methods
    % % % %     constructor method
    methods
        function obj = RGP(noInputs)
            %RGP Construct an instance of this class given no of inputs as
            %input argument
            obj.noInputs = noInputs;
        end
        
        % % % %         calculate covariance between 2 design points
        function val = calcCovariance(obj,x1,x2,hyperParams)
            % variance of latent function
            alphaSquared = hyperParams(1);
            % length scales
            lengthScales = hyperParams(2:end-1);
            % calculate covariance
            val = alphaSquared*exp(-0.5*(x1-x2)'*...
                (eye(length(lengthScales))./lengthScales.^2)*(x1-x2));
        end
        
        % % % %         calculate mean
        function val = calcMean(obj,x)
            % calculate mean
            val = 0;
        end
        
        % % % %         build covariance matrix and mean vector
        function [covMat,meanVec] = buildCovMatAndMeanVec(obj,x,hyperParams)
            % number of points = number of columns
            noPts = size(x,2);
            % initial matrices
            covMat = zeros(noPts);
            meanVec = NaN(noPts,1);
            % form lower triangular covariance matrix for covMat and
            % calculate mean
            for ii = 1:noPts
                for jj = ii:noPts
                    covMat(ii,jj) = obj.calcCovariance(x(:,ii),x(:,jj),...
                        hyperParams);
                end
                meanVec(ii,1) = obj.calcMean(x(:,ii));
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         calculate predicted mean and variance
        function [predMean,predVar] = calcPredMeanAndPredVar(obj,xt,x,y,...
                hyperParams)
            % number of trainings points
            noTP = size(x,2);
            % calculate mean vector and covariance matrix
            [covMat,meanVec] = obj.buildCovMatAndMeanVec(x,hyperParams);
            % add noise variance to Kmat
            Kx = covMat + hyperParams(end)*eye(noTP);
            % initialize vectors
            kStar = NaN(noTP,1);
            % populate kStar
            for ii = 1:noTP
                kStar(ii,1) = obj.calcCovariance(x(:,ii),xt,hyperParams);
            end
            % calculate kStarStar
            kStarStar = obj.calcCovariance(xt,xt,hyperParams);
            % calculate mStar
            mStar = obj.calcMean(xt);
            % calculate inverse of covariance matrix
            invK = inv(Kx);
            % predicted mean
            predMean = mStar + kStar'*invK*(y-meanVec);
            % predicted variance
            predVar = kStarStar - kStar'*invK*kStar;
            
        end
        
        % % % %         RGP regression
        function [muGt,cGt] = rgpRegression(obj,xt,yt,hyperParams,...
                initVals,muGt_1,cGt_1)
            % extract values from initialization structure
            x = initVals.basisVec;
            mX = initVals.initMeanVec;
            invkXX = initVals.initInvCovMat;
            % number of design/training points
            noTP = size(x,2);
            % calculate mean and covariance at candidate point
            mXt = obj.calcMean(xt);
            kXtXt = obj.calcCovariance(xt,xt,hyperParams);
            % calculate covariance of candidate wrt design points
            kXtX = NaN(1,noTP);
            for ii = 1:noTP
                kXtX(1,ii) = obj.calcCovariance(xt,x(:,ii),hyperParams);
            end
            % calculate Jt as per Huber Eqn. (8)
            Jt = kXtX*invkXX;
            % calculate B as per Huber Eqn. (7)
            B = kXtXt - Jt*kXtX';
            % calculate muP as per Huber Eqn. (6)
            muP = mXt + Jt*(muGt_1 - mX);
            % calculate cP as per Huber Eqn. (9)
            cP = B + Jt*cGt_1*Jt';
            % calculate Gt (kalman gain matrix) as per Huber Eqn. (12)
            Gt = cGt_1*Jt'*(cP + hyperParams(end))^-1;
            % calculate mean at step t as per Huber Eqn. (10)
            muGt = muGt_1 + Gt*(yt - muP);
            % calculate covariance at step t as per Huber Eqn. (11)
            cGt = cGt_1 - Gt*Jt*cGt_1;
        end
        
    end
end


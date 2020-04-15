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
            lengthScales = hyperParams(3:obj.noInputs+2);
            % calculate covariance
            val = alphaSquared*exp(-0.5*(x1-x2)'*...
                (eye(length(lengthScales))./lengthScales.^2)*(x1-x2));
        end
        
        % % % %         calculate mean
        function val = calcMean(obj,x)
            % calculate mean
            val = 0*(x'*x);
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
        
        % % % %         calculate marginal likelihood
        function logP = calcMarginalLikelihood(obj,x,y,hyperParams)
            % determine number of training points
            noTP = size(x,2);
            % build the covariance matrix
            kX = obj.buildCovMatAndMeanVec(x,hyperParams);
            % add signal noise to the covariance matrix
            kX = kX + eye(noTP)*hyperParams(2);
            % the data fit part
            dataFit = 0.5*y'*(kX\y);
            % complexity penalty
            complexPen = 0.5*log(det(kX));
            % normalizing constant
            normConstant = 0.5*noTP*log(2*pi);
            % marginal likelihood
            logP = - dataFit - complexPen - normConstant;
            
        end
        
        % % % %         calculate predicted mean and variance
        function [predMean,predVar] = calcPredMeanAndPredVar(obj,xt,x,y,...
                hyperParams)
            % number of trainings points
            noTP = size(x,2);
            % calculate mean vector and covariance matrix
            [covMat,meanVec] = obj.buildCovMatAndMeanVec(x,hyperParams);
            % add noise variance to Kmat
            Kx = covMat + hyperParams(2)*eye(noTP);
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
            % predicted mean
            predMean = mStar + kStar'*(Kx\(y-meanVec));
            % predicted variance
            predVar = kStarStar - kStar'*(Kx\kStar);
            
        end
        
        % % % %         RGP regression
        function [muGt,cGt] = rgpRegression(obj,xt,yt,hyperParams,...
                basisVec,initMeanVec,initInvCovMat,muGt_1,cGt_1)
            % extract values from initialization structure
            x = basisVec;
            mX = initMeanVec;
            invkXX = initInvCovMat;
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
            Gt = cGt_1*Jt'*(cP + hyperParams(2))^-1;
            % calculate mean at step t as per Huber Eqn. (10)
            muGt = muGt_1 + Gt*(yt - muP);
            % calculate covariance at step t as per Huber Eqn. (11)
            cGt = cGt_1 - Gt*Jt*cGt_1;
        end
        
                % % % %         RGP regression
        function [muGt,cGt,xn] = timeDependentRgpRegression(obj,xt,yt,...
                hyperParams,x,muGt_1,cGt_1)
            
            % number of design/training points
            noTP = size(x,2);            
            % build mean vectir and covariance matrix
            [kXX,mX] = obj.buildCovMatAndMeanVec(x,hyperParams);
            % calculate mean and covariance at candidate point
            mXt = obj.calcMean(xt);
            kXtXt = obj.calcCovariance(xt,xt,hyperParams);
            % calculate covariance of candidate wrt design points
            kXtX = NaN(1,noTP);
            for ii = 1:noTP
                kXtX(1,ii) = obj.calcCovariance(xt,x(:,ii),hyperParams);
            end
            % calculate Jt as per Huber Eqn. (8)
            Jt = kXtX/(kXX + (1e-5)*eye(noTP));
            % calculate B as per Huber Eqn. (7)
            B = kXtXt - Jt*kXtX';
            % calculate muP as per Huber Eqn. (6)
            muP = mXt + Jt*(muGt_1 - mX);
            % calculate cP as per Huber Eqn. (9)
            cP = B + Jt*cGt_1*Jt';
            % calculate Gt (kalman gain matrix) as per Huber Eqn. (12)
            Gt = cGt_1*Jt'*(cP + hyperParams(2))^-1;
            % calculate mean at step t as per Huber Eqn. (10)
            muGt = muGt_1 + Gt*(yt - muP);
            % calculate covariance at step t as per Huber Eqn. (11)
            cGt = cGt_1 - Gt*Jt*cGt_1;
            % add new measurement to the basis vector
            lia = ismember(x(1:end-1,:)',xt(1:end-1)','rows');
            lia = find(lia);
            xn = x;
            [~,minIdx] = min(xn(end,lia));
            xn(:,lia(minIdx)) = xt;
            
        end

        
    end
end


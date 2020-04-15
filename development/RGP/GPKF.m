classdef GPKF
    %GPKF Summary of this class goes here
    %   Detailed explanation goes here
    
    %% properties
    properties
        noSpatialIps
    end
    
    %% methods
    methods
        
        % % % %         constructor
        function obj = GPKF(noSpatialIps)
            obj.noSpatialIps = noSpatialIps;
        end
        
        % % % %         spatial kernel: squared exponential
        function val = meanFunction(obj,x)
            % % zero mean function
            val = 0*(x'*x)*obj.noSpatialIps;
        end
        
        % % % %         spatial kernel: squared exponential
        function val = spatialKernel(obj,s1,s2,covAmp,lengthScales)
            % % covariance equations
            val = covAmp*exp(-0.5*(s1-s2)'*(eye(obj.noSpatialIps)...
                ./lengthScales.^2)*(s1-s2));
        end
        
        % % % %         form covariance matrix and mean vector
        function covMat = buildSpatialCovMat(obj,x,covAmp,lengthScales)
            % number of points = number of columns
            noPts = size(x,2);
            % preallocate matricx
            covMat = zeros(noPts);
            % form lower triangular covariance matrix for covMat
            for ii = 1:noPts
                for jj = ii:noPts
                    covMat(ii,jj) = obj.spatialKernel(x(:,ii),x(:,jj),...
                        covAmp,lengthScales);
                end
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         temporal kernel: exponential
        function val = temporalKernel(obj,t1,t2,timeScale)
            % % covariance equation
            val = 1*exp(-timeScale*abs(t1-t2));
        end
        
        % % % %         calculate covaraince as a product of the two covariances
        function val = calcTotCovariance(obj,x1,x2,hyperParams)
            % % covariance amplitude or variance of latent function
            covAmp = hyperParams(1);
            % % length scales for spatial covariance
            lenScale = hyperParams(2:obj.noSpatialIps+1);
            % % time scale
            timeScale = hyperParams(obj.noSpatialIps+2);
            % % k = k_s*k_t
            val = obj.spatialKernel(x1(1:end-1),x2(1:end-1),covAmp,lenScale)...
                *obj.temporalKernel(x1(end),x2(end),timeScale);
        end
        
        % % % %         form covariance matrix and mean vector
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
                    covMat(ii,jj) = obj.calcTotCovariance(x(:,ii),x(:,jj),...
                        hyperParams);
                end
                meanVec(ii,1) = obj.meanFunction(x(:,ii));
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
            kX = kX + eye(noTP)*hyperParams(obj.noSpatialIps + 3);
            % the data fit part
            dataFit = 0.5*y'*(kX\y);
            % complexity penalty
            complexPen = 0.5*log(det(kX));
            % normalizing constant
            normConstant = 0.5*noTP*log(2*pi);
            % marginal likelihood
            logP = - dataFit - complexPen - normConstant;
        end
        
        % % % %         GPKF initialization
        function val = gpkfInitialize(obj,timeScale,timeStep)
            % % calculate F,H,Q as per Carron Eqn. (14)
            F = exp(-timeScale*timeStep);
            H = sqrt(2*1*timeScale);
            G = 1;
            fun = @(tau,timeScale) exp(-2*timeScale*tau);
            Q = integral(@(tau) fun(tau,timeScale),0,timeScale);
            % solve the Lyapunov equation for X
            sigma0 = (-G^2)/2*F;
            % outputs
            val.F = F;
            val.H = H;
            val.G = G;
            val.Q = Q;
            val.sigm0 = sigma0;
        end
        
        % % % %         GPKF implmentation
        function [predMean,predCov,skp1_kp1,ckp1_kp1] = ...
                gpkfRecurssion(obj,xDomain,xMeasure,sk_k,ck_k,Mk,yk,...
                Ks_12,F,Q,H,noiseVar)
            % total number of points in the entire domain of interest
            xDomainNP = size(xDomain,2);
            % number of measurable points which is subset of xDomain
            xMeasureNP = size(xMeasure,2);
            % number of points visited at each step which is a subset of xMeasure
            MkNP = size(Mk,2);
            % A matrix as per Carron conf. paper Eqn. (12)
            Amat = eye(xMeasureNP)*F;
            % Q matrix as per Carron conf. paper Eqn. (12)
            Qmat = eye(xMeasureNP)*Q;
            % H matrix as per Carron conf. paper Eqn. (12)
            Hmat = eye(xMeasureNP)*H;
            % R matrix as per Carron conf. paper Eqn. (12)
            Rmat = eye(MkNP)*noiseVar;
            % indicator matrix to find which points are visited at each iteration
            Ik = zeros(MkNP,xMeasureNP);
            % find points from xMeasure visited at iteration k
            lia = ismember(xMeasure',Mk','rows');
            lia = find(lia); % covert to numerical array instead of logical
            % populate the Ik matrix
            for ii = 1:MkNP
                Ik(ii,lia(ii)) = 1;
            end
            % C matrix as per Carron conf. paper Eqn. (12)
            Cmat = Ik*Ks_12*Hmat;
            % Kalman filter equations as per Carron conf. paper Eqn. (6)
            skp1_k = Amat*sk_k; % Eqn. (6a)
            ckp1_k = Amat*ck_k*Amat' + Qmat; % Eqn. (6b)
            Lkp1 = ckp1_k*Cmat'/(Cmat*ckp1_k*Cmat' + Rmat); % Eqn (6e)
            skp1_kp1 = skp1_k + Lkp1*(yk - Cmat*skp1_k); % Eqn (6c)
            ckp1_kp1 = ckp1_k - Lkp1*Cmat*ckp1_k; % Eqn (6d)
            % predicted value of mean and covarinace as per Todescato journal
            % paper Eqns. (13) and (14)
            predMean = Ks_12*Hmat*skp1_kp1; % Eqn. (13)
            predCov = Ks_12*Hmat*ckp1_kp1*Hmat'*Ks_12;
            % extend prediction over the entire domain
            if xDomainNP - xMeasureNP == 0
            end
            
        end
        
        
    end
end


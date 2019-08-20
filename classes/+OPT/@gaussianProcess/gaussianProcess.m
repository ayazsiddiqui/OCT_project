classdef gaussianProcess < dynamicprops
    %@GAUSSIANPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        noInputs
        kernelName
    end
    
    
    methods
        % contructor
        
        %                 function obj = gaussianProcess(obj)
        %                     %@GAUSSIANPROCESS Construct an instance of this class
        %                     %   Detailed explanation goes here
        %                     obj.noInputs = [];
        %                     obj.kernel = [];
        %                     obj.hyperParameters = [];
        %
        %                 end
        
        %% setters
        % set number of design variables
        function obj = set.noInputs(obj,val)
            if numel(val)== 1
                obj.noInputs = val;
            else
                error('Specify scalar value for number of inputs')
            end
        end
        
        % set kernel
        function obj = set.kernelName(obj,val)
            if ischar(val)
                obj.kernelName = val;
            else
                error('Kernel name should be a character string')
            end
        end
        
        % get kernel
        function obj = getkernel(obj)
            
            obj.addprop('kernel');
            
            switch obj.kernelName
                case 'squaredExponential'
                    obj.kernel = OPT.squaredExponential;
            end
            
            obj.kernel.noInputs = obj.noInputs;
        end
        
        %% other methods
        % objective function
        function val = objectiveFunction(obj,designPts)
            X = designPts;
            if size(designPts,1) ~= obj.noInputs
                error('Number of inputs mismatch')
            end
            val = -((X(1,:).^2 + X(2,:).^2)./50) + 1;
            val = reshape(val,[],1);
        end
        
        % build covariance matrix baed on kernel
        function val = buildCovarianceMatrix(obj,dsgnSet1,dsgnSet2)
            val = obj.kernel.buildCovarianceMatrix(dsgnSet1,dsgnSet2);
            
        end
        
        % calculate log likelihood
        function val = calcLogLikelihood(obj,dsgnSet,varargin)
            
            switch obj.kernelName
                case 'squaredExponential'
                    p = inputParser;
                    addParameter(p,'covarianceAmp',0,@isnumeric);
                    addParameter(p,'noiseVariance',0,@isnumeric);
                    addParameter(p,'lengthScale',1,@isnumeric);
                    
                    parse(p,varargin{:})
                    
                    obj.kernel.covarianceAmp = p.Results.covarianceAmp;
                    obj.kernel.noiseVariance = p.Results.noiseVariance;
                    obj.kernel.lengthScale = p.Results.lengthScale;
            end
            
            Kmat = obj.buildCovarianceMatrix(dsgnSet,dsgnSet);
            
            y = obj.objectiveFunction(dsgnSet);
            
            val = -1*(-0.5*(y'/Kmat*y) - 0.5*log(det(Kmat)));
        end
        
        % optimize hyper parameters
        function val = optimizeHyperParameters(obj,dsgnSet,initialGuess)
            
            A = []; b = [];
            Aeq = []; beq = [];
            
            lb = [0,0,zeros(1,obj.noInputs)];
            ub = [1e3*[1,1],1e3*ones(1,obj.noInputs)];
            
            switch obj.kernelName
                case 'squaredExponential'
                    val = fmincon(@(hyper) ...
                        obj.calcLogLikelihood(dsgnSet,...
                        'covarianceAmp',hyper(1),'noiseVariance',hyper(2),...
                        'lengthScale',hyper(3:end)),...
                        initialGuess,A,b,Aeq,beq,lb,ub);
                    
            end
        end
        
        % calculate predictive mean and variance
        function [predMean,predVar] = calcPredictiveMeanAndVariance(obj,postDsgn,trainDsgn,trainCovMat,trainFval)
            
            CovMatVecTranspose = obj.buildCovarianceMatrix(trainDsgn,postDsgn);
           
            % mean and variance
            nPost = size(postDsgn,2);
            
            muD = NaN(nPost,1);
            Var = NaN(nPost,1);
            for ii = 1:nPost
                muD(ii,1) = (CovMatVecTranspose(:,ii)'/trainCovMat)*trainFval;
                Var(ii,1) = obj.buildCovarianceMatrix(postDsgn(:,ii),postDsgn(:,ii))...
                    - (CovMatVecTranspose(:,ii)'/trainCovMat)*CovMatVecTranspose(:,ii);
            end
            
            predMean = muD;
            predVar = Var;
            
        end
        
        % calculate expected improvement
        function val = calcExpectedImprovement(obj,postDsgn,trainDsgn,trainCovMat,trainFval)
            
           [predMean,predVar] = obj.calcPredictiveMeanAndVariance(postDsgn,trainDsgn,trainCovMat,trainFval);
           fBest = max(trainFval);
           
           % http://krasserm.github.io/2018/03/21/bayesian-optimization/
           stdDev = sqrt(predVar);
           pd = makedist('Normal','mu',predMean,'sigma',stdDev);
           gm = gmdistribution(predMean,stdDev);
           
           if stdDev > 0
               Z = (predMean - fBest)/stdDev;
               val = -1*((predMean - fBest)*cdf(pd,Z) + stdDev*pdf(gm,Z));
           else
               val = 0;
           end
           
        end
        
    end
    
end




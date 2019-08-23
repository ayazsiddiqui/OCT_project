classdef gaussianProcess
    %@GAUSSIANPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        noInputs
        kernelName
        acquisitionFunctionName
        kernel
        acquisitionFunction
    end
    
    
    methods
        
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
                error('kernelName should be a character string')
            end
        end
        
        % set acquisition function
        function obj = set.acquisitionFunctionName(obj,val)
            if ischar(val)
                obj.acquisitionFunctionName = val;
            else
                error('acquisitionFunctionName should be a character string')
            end
        end
        
        % build GP
        function obj = buildKernel(obj)
            obj.kernel = kernels.(obj.kernelName);
            obj.kernel.noInputs = obj.noInputs;
        end
        
        function obj = buildAcquitisionFn(obj)
            obj.acquisitionFunction = acquisitionFunctions.(obj.acquisitionFunctionName);
        end
        
        %% other methods
        % objective function
        function val = objectiveFunction(obj,designPts)
            X = designPts;
            if size(designPts,1) ~= obj.noInputs
                error('Number of inputs mismatch')
            end
            %             % % % Park example 1
            %             val = -((X(1,:).^2 + X(2,:).^2)./50) + 1;
            % % % Park example 2
            val = 0.5*exp(-0.5*(X(2,:)-2).^2 - 0.5*(X(1,:)-2).^2)...
                +0.5*exp(-0.5*(X(1,:)+2).^2 - 0.5*(X(2,:)+2).^2);
            % % % https://www.hindawi.com/journals/mpe/2013/948303/ example
            %             val = exp(-((X(1,:)-4).^2 + (X(2,:)-4).^2)) + ...
            %                 exp(-((X(1,:)+4).^2 + (X(2,:)-4).^2)) + ...
            %                 2.*exp(-(X(1,:).^2 + X(2,:).^2)) + ...
            %                 2.*exp(-(X(1,:).^2 + (X(2,:)+4).^2));
            
            val = reshape(val,[],1);
        end
        
        % build covariance matrix baed on kernel
        function val = buildCovarianceMatrix(obj,dsgnSet1,dsgnSet2)
            val = obj.kernel.buildCovarianceMatrix(dsgnSet1,dsgnSet2);
            
        end
        
        % calculate log likelihood
        function val = calcLogLikelihood(obj,dsgnSet,dsgnFval,varargin)
            
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
            
            y = dsgnFval;
            
            val = 1*(-0.5*(y'/Kmat*y) - 0.5*log(det(Kmat)));
        end
        
        % optimize hyper parameters
        function val = optimizeHyperParameters(obj,dsgnSet,dsgnFval,initialGuess)
            
            A = []; b = [];
            Aeq = []; beq = [];
            
            % bounds
            lb = [eps,1e-2*ones(1,obj.noInputs)];
            ub = [10,10*ones(1,obj.noInputs)];
            nonlcon = [];
            options  = optimoptions('fmincon','Display','off');
            switch obj.kernelName
                case 'squaredExponential'
                    val = fmincon(@(hyper) ...
                        -obj.calcLogLikelihood(dsgnSet,dsgnFval,...
                        'covarianceAmp',hyper(1),'noiseVariance',obj.kernel.noiseVariance,...
                        'lengthScale',hyper(2:end)),...
                        initialGuess,A,b,Aeq,beq,lb,ub,nonlcon,options);
                    
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
        
        % calculate acquisition function
        function val = calcAcquisitionFunction(obj,testDsgn,bestFval,trainDsgn,trainCovMat,trainFval)
            
            [predMean,predVar] = obj.calcPredictiveMeanAndVariance(testDsgn,trainDsgn,trainCovMat,trainFval);
            
            switch obj.acquisitionFunctionName
                case 'expectedImprovement'
                    val = obj.acquisitionFunction.calcAcquisitionFunctionVal(predMean,predVar,bestFval);
                    
                case 'upperConfidenceBound'
                    val = obj.acquisitionFunction.calcAcquisitionFunctionVal(predMean,predVar);
                    
            end
            
        end
        
        % set bounds on desgin
        function val = calDesignBounds(obj,dsgnPt,tau,designLimits)
            
            lb = dsgnPt - tau;
            ub = dsgnPt + tau;
            
            lowLim = designLimits(:,1);
            hiLim = designLimits(:,2);
            
            belowLow = lb<lowLim;
            aboveHi = ub>hiLim;
            
            lb(belowLow) = lowLim(belowLow);
            ub(aboveHi) = hiLim(aboveHi);
            
            val = [lb(:) ub(:)];
        end
        
        % maximize acquisition function
        function [val,aFmax] = maximizeAcquisitionFunction(obj,bestFval,trainDsgn,trainCovMat,trainFval,initialPt,bounds)
            
            A = []; b = [];
            Aeq = []; beq = [];
            
            lb = bounds(:,1)';
            ub = bounds(:,2)';
            nonlcon = [];
            options  = optimoptions('fmincon','Display','off');
            [val,aFmax] = fmincon(@(optDsgn) ...
                -obj.calcAcquisitionFunction(optDsgn,bestFval,trainDsgn,trainCovMat,trainFval),...
                initialPt,A,b,Aeq,beq,lb,ub,nonlcon,options);
        end
        
        %% bayesian ascent
        
        function op = bayesianAscent(obj,trainDsgns,trainFval,trainOpHyp,iniPt,designLimits,iniTau,gamma,beta,maxIter)
            
            noIter = 1;
            
            % preallocate matrices
            noTrainDsgn = numel(trainFval);
            nt = noTrainDsgn;
            testDsgns = NaN(obj.noInputs,nt+maxIter);
            testFval = NaN(nt + maxIter,1);
            testOpHyp = NaN(size(trainOpHyp,1),maxIter);
            finPts = NaN(obj.noInputs,maxIter);
            finFval = NaN(maxIter,1);
            tau = NaN(size(iniTau,1),maxIter);
            % intitialize
            testDsgns(:,1:nt) = trainDsgns;
            testFval(1:nt,1) = trainFval;
            testOpHyp(:,1) = trainOpHyp;
            
            while noIter <= maxIter
                if noIter == 1
                    testDsgns(:,nt+noIter) = iniPt;
                    testFval(nt+noIter,1) = obj.objectiveFunction(iniPt);
                    testOpHyp(:,noIter) = trainOpHyp;
                    finPts(:,noIter) = iniPt;
                    finFval(noIter,1) = testFval(nt+noIter,1);
                    tau(:,noIter) = iniTau;
                    
                else
                    testDsgns(:,nt+noIter) = optPt;
                    testFval(nt+noIter,1) = optFval;
                    OpHyp = obj.optimizeHyperParameters(testDsgns(:,nt+noIter),testFval(nt+noIter,1),testOpHyp(:,noIter-1));
                    testOpHyp(:,noIter) = OpHyp;
                    finPts(:,noIter) = optPt;
                    finFval(noIter,1) = optFval;
                    
                    if finFval(noIter,1)-finFval(noIter-1,1) >= gamma*(1/noIter)*(max(finFval(1:noIter-1,1))-finFval(1))
                        tau(:,noIter) = beta*tau(:,noIter-1);
                        
                    else
                        tau(:,noIter) = iniTau;
                        fprintf('Bounds reset to initial bounds at iteration %d\n',noIter)
                    end
                end
                
                % step 1: optimizie hyper parameters
                obj.kernel.covarianceAmp = testOpHyp(1,noIter);
                obj.kernel.lengthScale = testOpHyp(2:end,noIter);
                
                % step 2: construct GP model
                testCovMat = obj.buildCovarianceMatrix(testDsgns(:,1:nt+noIter),testDsgns(:,1:nt+noIter));
                
                % select next input
                xLims = obj.calDesignBounds(finPts(:,noIter),tau(:,noIter),designLimits);
                
                % maximize acquisition function
                [optPt,~] = obj.maximizeAcquisitionFunction(max(finFval),testDsgns(:,1:nt+noIter),...
                    testCovMat,testFval(1:nt+noIter,1),finPts(:,noIter),xLims);
                optFval = obj.objectiveFunction(optPt);
                
                % convergence check
                noIter = noIter + 1;
                
            end
            
            % remove nan elemets
            chNaN = isnan(testDsgns);testDsgns(chNaN) = [];
            chNaN = isnan(testFval);testFval(chNaN) = [];
            chNaN = isnan(testOpHyp);testOpHyp(chNaN) = [];
            chNaN = isnan(finPts);finPts(chNaN) = [];
            chNaN = isnan(finFval);finFval(chNaN) = [];
            chNaN = isnan(tau);tau(chNaN) = [];
            
            op.testDsgns = testDsgns;
            op.testFval = testFval;
            op.testOpHyp = testOpHyp;
            op.finPts = [finPts,optPt];
            op.finFval = [finFval;optFval];
            op.tau = tau;
            
            
        end
        
        
    end
    
end




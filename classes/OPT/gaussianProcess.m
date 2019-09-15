classdef gaussianProcess
    %@GAUSSIANPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        noInputs
        kernel
        acquisitionFunction
    end
    
    
    methods
        %% contructor
        function obj = gaussianProcess(noInputs,varargin)
            p = inputParser;
            
            addRequired(p,'noInputs');
            addParameter(p,'kernel','squaredExponential',@ischar);
            addParameter(p,'acquisitionFunction','expectedImprovement',@ischar);
            
            parse(p,noInputs,varargin{:});
            
            obj.noInputs = p.Results.noInputs;
            obj.kernel = kernels.(p.Results.kernel);
            obj.kernel.noInputs = p.Results.noInputs;
            obj.acquisitionFunction = acquisitionFunctions.(p.Results.acquisitionFunction);
        end
        
        %% setters
        % set number of design variables
        function obj = set.noInputs(obj,val)
            if numel(val)== 1
                obj.noInputs = val;
            else
                error('Specify scalar value for number of inputs')
            end
        end
        
        
        %% other methods
        
        % build covariance matrix baed on kernel
        function val = buildCovarianceMatrix(obj,dsgnSet1,dsgnSet2)
            val = obj.kernel.buildCovarianceMatrix(dsgnSet1,dsgnSet2);
            
        end
        
        % calculate log likelihood
        function val = calcLogLikelihood(obj,dsgnSet,dsgnFval,varargin)
            
            switch class(obj.kernel)
                case 'kernels.squaredExponential'
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
            lb = [eps,1e-1*ones(1,obj.noInputs)]';
            ub = [10,10*ones(1,obj.noInputs)]';
            nonlcon = [];
            options  = optimoptions('fmincon','Display','off');
            switch class(obj.kernel)
                case 'kernels.squaredExponential'
                    val = fmincon(@(hyper) ...
                        -obj.calcLogLikelihood(dsgnSet,dsgnFval,...
                        'covarianceAmp',hyper(1,1),'noiseVariance',obj.kernel.noiseVariance,...
                        'lengthScale',hyper(2:end,1)),...
                        initialGuess,A,b,Aeq,beq,lb,ub,nonlcon,options);
                    %
                    %                     val = particleSwarmOpt(@(hyper) ...
                    %                         obj.calcLogLikelihood(dsgnSet,dsgnFval,...
                    %                         'covarianceAmp',hyper(1,1),'noiseVariance',obj.kernel.noiseVariance,...
                    %                         'lengthScale',hyper(2:end,1)),initialGuess,lb,ub,...
                    %                         'swarmSize',25,'cognitiveLR',0.4,'socialLR',0.2,'maxIter',20);
                    
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
            
            switch class(obj.acquisitionFunction)
                case 'acquisitionFunctions.expectedImprovement'
                    val = obj.acquisitionFunction.calcAcquisitionFunctionVal(predMean,predVar,bestFval);
                    
                case 'acquisitionFunctions.upperConfidenceBound'
                    val = obj.acquisitionFunction.calcAcquisitionFunctionVal(predMean,predVar);
                    
            end
            
        end
        
        % set bounds on desgin
        function [lowBound,upBound] = calDesignBounds(obj,dsgnPt,tau,designLimits)
            
            lowBound = dsgnPt - tau;
            upBound = dsgnPt + tau;
            
            lowLim = repmat(designLimits(:,1),1,size(dsgnPt,2));
            hiLim = repmat(designLimits(:,2),1,size(dsgnPt,2));
            
            belowLow = lowBound<lowLim;
            aboveHi = upBound>hiLim;
            
            lowBound(belowLow) = lowLim(belowLow);
            upBound(aboveHi) = hiLim(aboveHi);
            
            
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
        function [op,obj] = bayesianAscent(obj,trainDsgns,trainFval,finDsgns,finFval,...
                opHyp,tau,designLimits,iniTau,gamma,beta,noIter)
            
            % intitialize
            testDsgns = [trainDsgns,finDsgns(:,1:noIter)];
            testFval = [trainFval(:); finFval(1:noIter)];
            
            if noIter == 1
                iniGuess = ones(1+obj.noInputs,1);
                tau(:,noIter) = iniTau;
                
            else
                iniGuess = opHyp(:,noIter-1);
                
                if finFval(noIter,1)-finFval(noIter-1,1) >= gamma*(1/noIter)*(max(testFval(1:noIter-1,1))-finFval(1))
                    tau(:,noIter) = beta*tau(:,noIter-1);
                    
                else
                    tau(:,noIter) = iniTau;
                    fprintf('Bounds reset to initial bounds at iteration %d\n',noIter)
                end
            end
            
            % step 1: optimizie hyper parameters
            testOpHyp(:,noIter) = obj.optimizeHyperParameters(testDsgns,testFval,iniGuess);
            obj.kernel.covarianceAmp = testOpHyp(1,noIter);
            obj.kernel.lengthScale = testOpHyp(2:end,noIter);
            
            % step 2: construct GP model
            testCovMat = obj.buildCovarianceMatrix(testDsgns,testDsgns);
            
            % step 3: calc design bounds
            [lb,ub] = obj.calDesignBounds(finDsgns(:,noIter),tau(:,noIter),designLimits);
            xLims = [lb ub];
            
            % step 4: maximize acquisition function
            [optPt,optAq] = obj.maximizeAcquisitionFunction(max(finFval(1:noIter)),testDsgns,...
                testCovMat,testFval,finDsgns(:,noIter),xLims);
            
            [predMean,predVar] = obj.calcPredictiveMeanAndVariance(optPt,testDsgns,testCovMat,testFval);
            
            op.testDsgns = testDsgns;
            op.testCovMat = testCovMat;
            op.testFval = testFval;
            op.testOpHyp = testOpHyp;
            op.optPt = optPt;
            op.optAq = optAq;
            op.tau = tau;
            op.predMean = predMean;
            op.predVar = predVar;
            
        end
        
        function [op,obj] = mpcBayesianAscent(obj,trainDsgns,trainFval,finDsgns,finFval,...
                opHyp,tau,designLimits,iniTau,gamma,beta,noIter,predHorizon,ctrlHorizon)
            
            % intitialize
            
            
            if noIter == 1
                testDsgns = [trainDsgns,finDsgns(:,1:noIter)];
                testFval = [trainFval(:); finFval(1:noIter)];
                iniGuess = ones(1+obj.noInputs,1);
                tau(:,noIter) = iniTau;
                
            else
                testDsgns = [trainDsgns,finDsgns(:,1:1+(noIter-1)*ctrlHorizon)];
                testFval = [trainFval(:); finFval(1:1+(noIter-1)*ctrlHorizon)];
                iniGuess = opHyp(:,noIter-1);
                
                if finFval(noIter,1)-finFval(noIter-1,1) >= gamma*(1/noIter)*(max(testFval(1:noIter-1,1))-finFval(1))
                    tau(:,noIter) = beta*tau(:,noIter-1);
                    
                else
                    tau(:,noIter) = iniTau;
                    fprintf('Bounds reset to initial bounds at iteration %d\n',noIter)
                end
            end
            
            % step 1: optimizie hyper parameters
            testOpHyp(:,noIter) = obj.optimizeHyperParameters(testDsgns,testFval,iniGuess);
            obj.kernel.covarianceAmp = testOpHyp(1,noIter);
            obj.kernel.lengthScale = testOpHyp(2:end,noIter);
            
            % step 2: construct GP model
            testCovMat = obj.buildCovarianceMatrix(testDsgns,testDsgns);
            
            % step 3: calc design bounds
            tauMpc = NaN(obj.noInputs,predHorizon);
            
            for ii = 1:predHorizon
                tauMpc(:,ii) = ii*tau(:,noIter);
            end
            
            % step 4: maximize acquisition function
            [lb,ub] = obj.calDesignBounds(repmat(finDsgns(:,1+(noIter-1)*ctrlHorizon),1,predHorizon),tauMpc,designLimits);
            
            iniGuess = repmat(finDsgns(:,1+(noIter-1)*ctrlHorizon),1,predHorizon) + 0.25*(ub-lb);
            
            
%                         [mpcOptPts,mpcOptFval] = particleSwarmOpt(@(pDsgn) obj.mpcPrediction...
%                             (pDsgn,testFval,testDsgns,testCovMat,testFval),iniGuess,lb,ub,...
%                             'swarmSize',1000,'cognitiveLR',0.4,'socialLR',0.2,'maxIter',15);
%             
%                         if isa(obj.acquisitionFunction,'acquisitionFunctions.expectedImprovement')
%                             thres = 1e-3;
%                             if mpcOptFval <thres
%                                 [~,idx] = max(testFval);
%                                 mpcOptPts = testDsgns(:,idx);
%                                 fprintf('Maximum expected improvement less than user defined threshold of %6.6f.\n',thres)
%             
%                             end
%                         end
%             
            [mpcOptPts,mpcOptFval] = sequentialParticleSwarmOpt(@(pDsgn) obj.mpcPrediction...
                (pDsgn,finFval,testDsgns,testCovMat,testFval),finDsgns(:,noIter),predHorizon,designLimits(:,1),designLimits(:,2),...
                'swarmSize',5000,'cognitiveLR',0.1,'socialLR',0.2,'maxIter',3,'stepPerc',0.2);
            
            % outputs
            op.mpcOptPts = mpcOptPts;
            op.mpcOptFval = mpcOptFval;
            
            [mpcPredMean,mpcPredVar] = obj.calcPredictiveMeanAndVariance(mpcOptPts,testDsgns,testCovMat,testFval);
            
            optPt = mpcOptPts(:,1:ctrlHorizon);
            optAq = mpcOptFval;
            
            op.testDsgns = testDsgns;
            op.testCovMat = testCovMat;
            op.testFval = testFval;
            op.testOpHyp = testOpHyp;
            op.optPt = optPt;
            op.optAq = optAq;
            op.tau = tau;
            op.mpcPredMean = mpcPredMean;
            op.mpcPredVar = mpcPredVar;
            op.mpcOptPts = mpcOptPts;
            op.mpcOptFval = mpcOptFval;
            
        end
        
        function val = mpcPrediction(obj,postDsgns,finFval,testDsgns,testCovMat,testFval)
            
            AqVals = obj.calcAcquisitionFunction(postDsgns,max(finFval),testDsgns,testCovMat,testFval);
            nElem = numel(AqVals);
%             val = sum(AqVals(:));
            %             val = AqVals(end);
            val = sum([1:nElem]'.*AqVals);
            %             val = sum([nElem:-1:1]'.*AqVals);
            
            
        end
        
    end
    
end




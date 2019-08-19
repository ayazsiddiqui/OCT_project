classdef gaussianProcess
    %@GAUSSIANPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        noInputs
        kernel
        hyperParameters
    end
    
    methods
        %% contructor
        
        %         function obj = gaussianProcess(inputArg1,inputArg2)
        %             %@GAUSSIANPROCESS Construct an instance of this class
        %             %   Detailed explanation goes here
        %             obj.Property1 = inputArg1 + inputArg2;
        %         end
        
        %% setters
        % set number of design variables
        function obj = setNoInputs(obj,val)
            obj.noInputs = val;
        end
        
        % set kernel
        function obj = setKernel(obj,val)
            if ischar(val)
                obj.kernel = val;
            else
                error('Kernel should be a character string')
            end
        end
        
        % set hyperParameters
        function obj = setHyperParameters(obj,val)
            obj.hyperParameters = val;
        end
        
        %% other methods
        % objective function (Ideally would be an input)
        function val = objectiveFunction(obj,designPts)
            X = designPts;
            if size(designPts,1) ~= obj.noInputs
                error('Number of inputs mismatch')
            end
            val = -((X(1,:).^2 + X(2,:).^2)./50) + 1;
            val = reshape(val,[],1);
        end
        
        % build covariance mat
        function covMat = buildCovarianceMatrix(obj,dsgnSet1,dsgnSet2,varargin)
            
            % function that builds covariance matrix between design sets
            % inputs datapoints, covariance amplitude, noise variance and length scale
            
            % length scale quantifies the relevance of the components in the input
            % vector for predicting the outputs
            
            switch obj.kernel
                case 'squaredExponential'
                    % squared exponential
                    cov = @(x1,x2,covAmp,noiseVar,lengthScl)...
                        covAmp*(exp(-0.5*((x1-x2)'*(eye(numel(lengthScl))./(lengthScl.^2))*(x1-x2))));
            end
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
            
            if isequal(dsgnSet1,dsgnSet2)
                for ii = 1:covMatSize(1)
                    for jj = ii:covMatSize(2)
                        covMat(ii,jj) = cov(dsgnSet1(:,ii),dsgnSet2(:,jj),...
                            p.Results.covAmplitude,p.Results.noiseVariance,p.Results.lengthScale);
                        
                    end
                end
                covMat = covMat + covMat' - eye(length(covMat)).*diag(covMat - p.Results.noiseVariance);
                
            else
                for ii = 1:covMatSize(1)
                    for jj = 1:covMatSize(2)
                        covMat(ii,jj) = cov(dsgnSet1(:,ii),dsgnSet2(:,jj),...
                            p.Results.covAmplitude,p.Results.noiseVariance,p.Results.lengthScale);
                        %             if ii == jj
                        %                 covMat(ii,jj) = covMat(ii,jj) + p.Results.noiseVariance;
                        %             end
                        
                    end
                end
            end
            
            
            
        end
    end
    
end


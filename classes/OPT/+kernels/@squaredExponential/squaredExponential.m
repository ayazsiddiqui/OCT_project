classdef squaredExponential
    %SQUAREDEXPONENTIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        noInputs
        covarianceAmp
        noiseVariance
        lengthScale
    end
    
    methods
        %% setter
        % set number of inputs
        function obj = set.noInputs(obj,val)
            if numel(val) == 1
                obj.noInputs = val;
            else
                error('Specify scalar value for number of inputs')
            end
        end
        
        % set covariance amplitutde
        function obj = set.covarianceAmp(obj,val)
            if numel(val) == 1
                obj.covarianceAmp = val;
            else
                error('Specify scalar value for covariance amplitude');
            end
        end
        
        % set noise covariance
        function obj = set.noiseVariance(obj,val)
            if numel(val) == 1
                obj.noiseVariance = val;
            else
                error('Specify scalar value for noise covariance');
            end
        end
        
        function obj = set.lengthScale(obj,val)
            if isequal(size(val),[obj.noInputs,1])
                obj.lengthScale = val;
            else
                error('Incorrect size of length scale matrix')
            end
        end
        
        %% other methods
        % build covariance mat
        function val = buildCovarianceMatrix(obj,dsgnSet1,dsgnSet2)
            
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
                            obj.covarianceAmp,obj.noiseVariance,obj.lengthScale);
                        
                    end
                end
                covMat = covMat + covMat' - eye(length(covMat)).*diag(covMat - obj.noiseVariance);
                
            else
                for ii = 1:covMatSize(1)
                    for jj = 1:covMatSize(2)
                        covMat(ii,jj) = cov(dsgnSet1(:,ii),dsgnSet2(:,jj),...
                            obj.covarianceAmp,obj.noiseVariance,obj.lengthScale);
                    end
                end
            end
            
            val = covMat;
            
        end
        
        
    end
end


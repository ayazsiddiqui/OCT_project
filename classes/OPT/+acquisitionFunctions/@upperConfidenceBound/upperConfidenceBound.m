classdef upperConfidenceBound
    %UPPERCONFIDENCEBOUND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        explorationFactor
    end
    
    % methods
    methods
        
        %% set mehods
        function obj = set.explorationFactor(obj,val)
            if numel(val) == 1 && val>= 0
                obj.explorationFactor = val;
            else
                error('Bad value for explorationFactor');
            end
        end
        
        %% other methods
        function val = calcAcquisitionFunctionVal(obj,predMean,predVar)
            val = 1*(predMean + obj.explorationFactor*predVar);
        end
    end
end


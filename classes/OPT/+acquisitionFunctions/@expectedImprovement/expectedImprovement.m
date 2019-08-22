classdef expectedImprovement
    %EXPECTEDIMPROVEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        
        function val = calcAcquisitionFunctionVal(obj,predMean,predVar,bestFval)
            % http://krasserm.github.io/2018/03/21/bayesian-optimization/
            fBest = bestFval;
            
            stdDev = sqrt(predVar);
            pd = makedist('Normal','mu',0,'sigma',1);
            gm = gmdistribution(0,1);
            
            if stdDev>0
                Z = (predMean-fBest)/stdDev;
                val = 1*(((predMean-fBest)*cdf(pd,Z)) + stdDev*pdf(gm,Z));
            else
                val = 0;
            end
        end
    end
end


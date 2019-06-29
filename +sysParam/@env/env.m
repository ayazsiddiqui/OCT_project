classdef env
    %ENV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gravAccel
        flowDensity
        inertialFlowVel
    end
    
    %% constructor
    methods
        function thisEnv = env()
            thisEnv.gravAccel.value = 9.81;
            thisEnv.flowDensity.value = 1000;
            thisEnv.inertialFlowVel.value = [];
        end
        
    end
    
    
end


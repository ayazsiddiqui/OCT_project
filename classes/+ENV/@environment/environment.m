classdef environment
    %ENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        fluidDensity
        inertialFlowVel
    end
    
    properties (Constant)
        gravAccel = SIM.parameter('Value',9.81,'Unit','m/s^2','Description','Graviational Acceleration');
    end
    
    
    methods
        %% constructor
        function obj = environment
            %ENVIRONMENT Construct an instance of this class
            obj.fluidDensity = SIM.parameter('Unit','kg/m^3','Description','Fluid density');
            obj.inertialFlowVel = SIM.parameter('Unit','m/s','Description','Inertial constant flow velocity');
        end
        
        %% setters
        function setFluidDensity(obj,val,units)
            obj.fluidDensity.setValue(val,units);
        end
        
        function setInertialFlowVel(obj,val,units)
            obj.inertialFlowVel.setValue(reshape(val,3,1),units);
        end
        
        
        %% other methods
        
        % scale environment
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
        
    end
end


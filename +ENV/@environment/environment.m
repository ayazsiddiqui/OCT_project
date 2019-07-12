classdef environment
    %ENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        inertialFlowVel
    end
    
    properties (Constant)
       gravAccel = SIM.parameter('Value',9.81,'Unit','m/s^2','Description','Graviational Acceleration');
    end
    
    properties (Dependent)
        fluidDensity
    end
    
    methods
        %% constructor
        function obj = environment
            %ENVIRONMENT Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.inertialFlowVel = SIM.parameter('Unit','m/s','Description','Inertial constant flow velocity');
        end
        
        %% setters
        function setLengthScale(obj,val,units)
            obj.lengthScale.setValue(val,units);
        end
        
        function setDensityScale(obj,val,units)
            obj.densityScale.setValue(val,units);
        end
        
        function setInertialFlowVel(obj,val,units)
           obj.inertialFlowVel.setValue(val,units);
        end
        
        %% getters
        function val = get.fluidDensity(obj)
            val =  SIM.parameter('Value',1000*obj.densityScale.Value,...
                'Unit','kg/m^3','Description','Fluid density');
        end
        
        %% other methods
        
        % scale environment
        function scaleEnvironment(obj)
            LS = obj.lengthScale.Value;
            
            obj.setInertialFlowVel(obj.inertialFlowVel.Value.*LS^0.5,'m/s');
        end

        
        
    end
end


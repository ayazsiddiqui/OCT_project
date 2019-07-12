classdef turbine
    %TURBINE Summary of this class goes here
    
    properties
        lengthScale
        densityScale
        numTurbines
        turbDiameter
        turbDragCoeff
        turbPowerCoeff
    end
    
    methods
        %% constructor
        function obj = turbine
            %TURBINE Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTurbines = SIM.parameter('Description','Number of turbines');
            obj.turbDiameter = SIM.parameter('Unit','m','Description','Turbine Diameter');
            obj.turbDragCoeff = SIM.parameter('Description','Turbine drag coefficient');
            obj.turbPowerCoeff = SIM.parameter('Description','Turbine Power Coeff');
            
        end
        
        %% setters
        function setLengthScale(obj,val,units)
            obj.lengthScale.setValue(val,units);
        end
        
        function setDensityScale(obj,val,units)
            obj.densityScale.setValue(val,units);
        end
        
        function setNumTurbines(obj,val,units)
            obj.numTurbines.setValue(val,units);
        end
        
        function setTurbDiameter(obj,val,units)
            val = reshape(val,1,[]);
            obj.turbDiameter.setValue(val,units)
           if length(val) ~= obj.numTurbines.Value
               error('Number of given values not equal to number of turbines')
           end
               
        end
        
        function setTurbDragCoeff(obj,val,units)
            val = reshape(val,1,[]);
            obj.turbDragCoeff.setValue(val,units)
           if length(val) ~= obj.numTurbines.Value
               error('Number of given values not equal to number of turbines')
           end
               
        end
        
        function setTurbPowerCoeff(obj,val,units)
            val = reshape(val,1,[]);
            obj.turbPowerCoeff.setValue(val,units)
            if length(val) ~= obj.numTurbines.Value
                error('Number of given values not equal to number of turbines')
            end
            
        end
        
        %% other methods
        
        % scale turbine
        function scaleTurbine(obj)
            LS = obj.lengthScale.Value;

            % scale turbine diameter
            obj.setTurbDiameter(obj.turbDiameter.Value*LS,'m');
            
        end
        
    end
end


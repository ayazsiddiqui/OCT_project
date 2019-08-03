classdef turbine
    %TURBINE Summary of this class goes here
    
    properties (SetAccess = private)
        numTurbines
        turbDiameter
        turbDragCoeff
        turbPowerCoeff
    end
    
    methods
        %% constructor
        function obj = turbine
            %TURBINE Construct an instance of this class
            obj.numTurbines = SIM.parameter('Description','Number of turbines');
            obj.turbDiameter = SIM.parameter('Unit','m','Description','Turbine Diameter');
            obj.turbDragCoeff = SIM.parameter('Description','Turbine drag coefficient');
            obj.turbPowerCoeff = SIM.parameter('Description','Turbine Power Coeff');
            
        end
        
        %% setters
        function setNumTurbines(obj,val,units)
            obj.numTurbines.setValue(val,units);
        end
        
        function setTurbDiameter(obj,val,units)
           if numel(val) ~= obj.numTurbines.Value
               error('Number of given values not equal to number of turbines');
           else
               obj.turbDiameter.setValue(reshape(val,1,[]),units);
           end
               
        end
        
        function setTurbDragCoeff(obj,val,units)
           if numel(val) ~= obj.numTurbines.Value
               error('Number of given values not equal to number of turbines');
           else
               obj.turbDragCoeff.setValue(reshape(val,1,[]),units);
           end
               
        end
        
        function setTurbPowerCoeff(obj,val,units)
            if numel(val) ~= obj.numTurbines.Value
                error('Number of given values not equal to number of turbines');
            else
                obj.turbPowerCoeff.setValue(reshape(val,1,[]),units);
            end
            
        end
        
        %% other methods
        
        % scale turbine
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
    end
end


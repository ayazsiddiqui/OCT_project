classdef winch
    %WINCH Summary of this class goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        wnchMaxTugSpeed
        wnchMaxReleaseSpeed
        wnchTimeConstant
    end
    
    methods
        %% constructor
        function obj = winch
            %WINCH Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTethers  = SIM.parameter('Description','Number of tethers');
            obj.wnchMaxTugSpeed = SIM.parameter('Unit','m/s','Description','Max winch pull in (wind up) speed');
            obj.wnchMaxReleaseSpeed = SIM.parameter('Unit','m/s','Description','Max winch push out (wind out) speed');
            obj.wnchTimeConstant = SIM.parameter('Unit','s','Description','Winch time constant');
            
        end
        
        %% setters
        function setLengthScale(obj,val,units)
            obj.lengthScale.setValue(val,units);
        end
        
        function setDensityScale(obj,val,units)
            obj.densityScale.setValue(val,units);
        end
        
        function setNumTethers(obj,val,units)
            obj.numTethers.setValue(val,units);
        end
        
        function setWnchMaxTugSpeed(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.wnchMaxTugSpeed.setValue(-abs(val),units);
            end
        end
        
        function setWnchMaxReleaseSpeed(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.wnchMaxReleaseSpeed.setValue(val,units);
            end
        end
        
        function setWnchTimeConstant(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.wnchTimeConstant.setValue(val,units);
            end
        end
        
        %% other methods
        
        % scale winches
        function scaleWinch(obj)
           LS = obj.lengthScale.Value;
           
           obj.setWnchMaxTugSpeed(obj.wnchMaxTugSpeed.Value.*LS^0.5,'m/s');
           obj.setWnchMaxReleaseSpeed(obj.wnchMaxReleaseSpeed.Value.*LS^0.5,'m/s');
           obj.setWnchTimeConstant(obj.wnchTimeConstant.Value.*LS^0.5,'s');
           
        end
            
        
    end
end


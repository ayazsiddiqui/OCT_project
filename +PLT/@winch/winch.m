classdef winch
    %WINCH Summary of this class goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        wnchMaxTugSpeed
        wnchMaxReleaseSpeed
        wnchTimeConstant
        initThrLength
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
            obj.initThrLength = SIM.parameter('Unit','m','Description','Initial tether length');
            
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
                obj.wnchMaxTugSpeed.setValue(-abs(reshape(val,1,[])),units);
            end
        end
        
        function setWnchMaxReleaseSpeed(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.wnchMaxReleaseSpeed.setValue(reshape(val,1,[]),units);
            end
        end
        
        function setWnchTimeConstant(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.wnchTimeConstant.setValue(reshape(val,1,[]),units);
            end
        end
        
        function setInitThrLength(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.initThrLength.setValue(reshape(val,1,[]),units);
            end
        end
        
        %% other methods
        
        % scale winches
        function scaleWinch(obj)
           LS = obj.lengthScale.Value;
           
           obj.setWnchMaxTugSpeed(obj.wnchMaxTugSpeed.Value.*LS^0.5,'m/s');
           obj.setWnchMaxReleaseSpeed(obj.wnchMaxReleaseSpeed.Value.*LS^0.5,'m/s');
           obj.setWnchTimeConstant(obj.wnchTimeConstant.Value.*LS^0.5,'s');
%            obj.setInitThrLength(obj.initThrLength.Value.*LS,'m');
           
        end
        
        
        % initial tether length
        function calcInitTetherLength(obj,vehicle,gndStn,environment)
            % dummy vars
            nt = obj.numTethers.Value;
            Rn_1 = NaN(3,nt);
            init_L = NaN(1,nt);
            
            for ii = 1:nt
                Rn_1(:,ii) = vehicle.init_inertialCmPos.Value + ...
                    rotation_sequence(vehicle.init_euler.Value)*...
                    vehicle.thrAttchPts.Value(:,ii) - ...
                    ( rotation_sequence([0;0;gndStn.init_euler.Value])*...
                    gndStn.thrAttchPts.Value(:,ii) );
                init_L(1,ii) = norm(Rn_1(:,ii));
            end
            
            obj.setInitThrLength(init_L,'m');
            
            
        end
            
            
        
            
        
    end
end


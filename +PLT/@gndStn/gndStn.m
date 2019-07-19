classdef gndStn
    %GNDSTN Summary of this class goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        Izz
        dampingCoeff
        % freeSpeen
        freeSpinSwitch
        % tether attachment points
        thrAttchPts
        % initial conditions
        init_euler
        init_angVel
    end
    
    methods
        %% constructor
        function obj = gndStn
            %GNDSTN Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTethers   = SIM.parameter('Description','Number of tethers');
            obj.Izz          = SIM.parameter('Unit','kg*m^2','Description','Izz');
            obj.dampingCoeff = SIM.parameter('Unit','N*m*s','Description','Ground station damping coeff');
            obj.freeSpinSwitch = SIM.parameter('Value',0,'Unit','','Description',...
                'Switch to enable or disable free spinning: 1=ON, 0=OFF');
            obj.thrAttchPts  = SIM.parameter('Unit','m','Description','Ground station tether attachment points');
            % initial conditions
            obj.init_euler         = SIM.parameter('Value',0,'Unit','rad','Description','Initial Euler angles');
            obj.init_angVel        = SIM.parameter('Value',0,'Unit','rad/s','Description','Initial angular velocities');
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
        
        function setIzz(obj,val,units)
            obj.Izz.setValue(val,units);
        end
        
        function setDampingCoeff(obj,val,units)
            obj.dampingCoeff.setValue(val,units);
        end
        
        function setFreeSpinSwitch(obj,val,units)
            if val~=1 && val~=0
                error('Invalid value, enter 0 or 1');
            else
                obj.freeSpinSwitch.setValue(val,units);
            end
        end
        
        % tether attachment points
        function setThrAttchPts(obj,vehicle)
                obj.thrAttchPts.setValue(vehicle.thrAttchPts.Value.*...
                    [ones(2,size(vehicle.thrAttchPts.Value,2));...
                zeros(1,size(vehicle.thrAttchPts.Value,2))],'m');
        end
        
        % initial conditions     
        function setInitialEuler(obj,val,units)
            if length(val) ~= 1
                error('Invalid matrix size, enter a single scalar value')
            end
            obj.init_euler.setValue(val,units);
        end
        
        function setInitialAngVel(obj,val,units)
            if length(val) ~= 1
                error('Invalid matrix size, enter a single scalar value')
            end
            obj.init_angVel.setValue(val,units);
        end
        
        %% other methods
        function scaleGndStn(obj)
            LS = obj.lengthScale.Value;
            DS = obj.densityScale.Value;
            
            % scale inertias and damping coeff
            obj.setIzz(obj.Izz.Value*(DS*LS^5),'kg*m^2');
            obj.setDampingCoeff(obj.dampingCoeff.Value*(DS*LS^4.5),'N*m*s');
            
            % scale initial conditions
            obj.setInitialAngVel(obj.init_angVel.Value.*(1/LS^0.5),'rad/s');

        end
        
        
    end
end


classdef gndStn
    %GNDSTN Summary of this class goes here
    
    properties (SetAccess = private)
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
            obj.numTethers   = SIM.parameter('Description','Number of tethers');
            obj.Izz          = SIM.parameter('Unit','kg*m^2','Description','Izz');
            obj.dampingCoeff = SIM.parameter('Unit','N*m*s','Description','Ground station damping coeff');
            obj.freeSpinSwitch = SIM.parameter('Value',0,'Unit','','Description',...
                'Switch to enable or disable free spinning: 1=ON, 0=OFF');
            obj.thrAttchPts  = SIM.parameter('Unit','m','Description','Ground station tether attachment points');
            % initial conditions
            obj.init_euler         = SIM.parameter('Unit','rad','Description','Initial Euler angles');
            obj.init_angVel        = SIM.parameter('Unit','rad/s','Description','Initial angular velocities');
        end
        
        %% setters
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
        function setThrAttchPts(obj,val,units)
            if numel(val) ~= 3*obj.numTethers.Value
                error('Number of attachment point vectors not equal to number of tethers');
            else
                obj.thrAttchPts.setValue(val,units);
            end
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
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
        
    end
end


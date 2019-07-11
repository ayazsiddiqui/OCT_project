classdef gndStn
    %GNDSTN Summary of this class goes here
    
    properties
        lengthScale
        densityScale
        numTethers
        Izz
        dampingCoeff
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
        
        % tether attachment points
        function setThrAttchPts(obj,val,units)
            obj.thrAttchPts.setValue(val,units)
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
            
            % scale inertias and damping coeff
            obj.setIzz(obj.Izz.Value*LS^5,'kg*m^2');
            obj.setDampingCoeff(obj.dampingCoeff.Value*LS^4.5,'N*m*s');
            
            % scale tether attachment vector
            obj.setThrAttchPts(obj.thrAttchPts.Value.*LS,'m');

            % scale initial conditions
            obj.setInitialAngVel(obj.init_angVel.Value.*(1/LS^0.5),'rad/s');

        end
        
        
    end
end


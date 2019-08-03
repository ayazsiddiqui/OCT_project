classdef oneThrCtlr
    %ONETHRCTLR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        % aileron control
        aileronKp
        aileronKi
        aileronKd
        aileronTau
        aileronMaxDef
        % elevator control
        elevatorKp
        elevatorKi
        elevatorKd
        elevatorTau
        elevatorMaxDef
        % rudder control
        rudderKp
        rudderKi
        rudderKd
        rudderTau
        rudderMaxDef
    end
    
    properties (Constant)
        CsAllocationMat = SIM.parameter('Value',[-1 0 0; 1 0 0; 0 -1 0; 0 0 1],...
            'Unit','','Description','Matrix that allocates control surface controller errors to CS deflections');
    end
    
    methods
        %% constructor
        function obj = oneThrCtlr
            %ONETHRCTLR Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTethers   = SIM.parameter('Description','Number of tethers');
            % aileron control
            obj.aileronKp      = SIM.parameter('Unit','deg/deg','Description','Aileron Kp');
            obj.aileronKi      = SIM.parameter('Unit','deg/(deg*s)','Description','Aileron Ki');
            obj.aileronKd      = SIM.parameter('Unit','deg/(deg/s)','Description','Aileron Kd');
            obj.aileronTau     = SIM.parameter('Unit','s','Description','Aileron time constant');
            obj.aileronMaxDef  = SIM.parameter('Unit','deg','Description','Aileron max deflection');
            % elevator control
            obj.elevatorKp     = SIM.parameter('Unit','deg/deg','Description','Elevator Kp');
            obj.elevatorKi     = SIM.parameter('Unit','deg/(deg*s)','Description','Elevator Ki');
            obj.elevatorKd     = SIM.parameter('Unit','deg/(deg/s)','Description','Elevator Kd');
            obj.elevatorTau    = SIM.parameter('Unit','s','Description','Elevator time constant');
            obj.elevatorMaxDef = SIM.parameter('Unit','deg','Description','Elevator max deflection');
            % rudder control
            obj.rudderKp       = SIM.parameter('Unit','deg/deg','Description','Rudder Kp');
            obj.rudderKi       = SIM.parameter('Unit','deg/(deg*s)','Description','Rudder Kp');
            obj.rudderKd       = SIM.parameter('Unit','deg/(deg/s)','Description','Rudder Kp');
            obj.rudderTau      = SIM.parameter('Unit','s','Description','Rudder time constant');
            obj.rudderMaxDef   = SIM.parameter('Unit','deg','Description','Rudder max deflection');
            
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
        
        % set aileron gains
        function setAileronKp(obj,val,units)
            obj.aileronKp.setValue(val,units);
        end
        
        function setAileronKi(obj,val,units)
            obj.aileronKi.setValue(val,units);
        end

        function setAileronKd(obj,val,units)
            obj.aileronKd.setValue(val,units);
        end
        
        function setAileronTau(obj,val,units)
            obj.aileronTau.setValue(val,units);
        end
        
        function setAileronMaxDef(obj,val,units)
           obj.aileronMaxDef.setValue(val,units) 
        end
        
        % set elevator gains
        function setElevatorKp(obj,val,units)
            obj.elevatorKp.setValue(val,units);
        end
        
        function setElevatorKi(obj,val,units)
            obj.elevatorKi.setValue(val,units);
        end

        function setElevatorKd(obj,val,units)
            obj.elevatorKd.setValue(val,units);
        end
        
        function setElevatorTau(obj,val,units)
            obj.elevatorTau.setValue(val,units);
        end
        
        function setElevatorMaxDef(obj,val,units)
           obj.elevatorMaxDef.setValue(val,units) 
        end
        
        % set rudder gains
        function setRudderKp(obj,val,units)
            obj.rudderKp.setValue(val,units);
        end
        
        function setRudderKi(obj,val,units)
            obj.rudderKi.setValue(val,units);
        end

        function setRudderKd(obj,val,units)
            obj.rudderKd.setValue(val,units);
        end
        
        function setRudderTau(obj,val,units)
            obj.rudderTau.setValue(val,units);
        end
        
        function setRudderMaxDef(obj,val,units)
           obj.rudderMaxDef.setValue(val,units) 
        end
        
        %% other methods
        
        % scale controller
        function scaleOneThrCtlr(obj)
            
            LS = obj.lengthScale.Value;
            % scale aileron gains
            obj.setAileronKi(obj.aileronKi.Value*(1/LS^0.5),'deg/(deg*s)');
            obj.setAileronKd(obj.aileronKd.Value*LS^0.5,'deg/(deg/s)');
            obj.setAileronTau(obj.aileronTau.Value*LS^0.5,'s');
            
            % scale elevator gains
            obj.setElevatorKi(obj.elevatorKi.Value*(1/LS^0.5),'deg/(deg*s)');
            obj.setElevatorKd(obj.elevatorKd.Value*LS^0.5,'deg/(deg/s)');
            obj.setElevatorTau(obj.elevatorTau.Value*LS^0.5,'s');
            
            % scale rudder gains
            obj.setRudderKi(obj.rudderKi.Value*(1/LS^0.5),'deg/(deg*s)');
            obj.setRudderKd(obj.rudderKd.Value*LS^0.5,'deg/(deg/s)');
            obj.setRudderTau(obj.rudderTau.Value*LS^0.5,'s');
            
        end


    end
end


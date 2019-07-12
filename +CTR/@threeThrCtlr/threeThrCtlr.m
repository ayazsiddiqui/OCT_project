classdef threeThrCtlr
    %THREETHRCTLR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        % altitude tether control
        altiTetherKp
        altiTetherKi
        altiTetherKd
        altiTetherTau
        altiErrorSat
        % pitch tether control
        pitchTetherKp
        pitchTetherKi
        pitchTetherKd
        pitchTetherTau
        % roll tether control
        rollTetherKp
        rollTetherKi
        rollTetherKd
        rollTetherTau
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
       thrAllocationMat = SIM.parameter('Value',[1 .5 -.5; 1 -.5 0; 1 .5 .5],...
           'Unit','','Description','Matrix that allocates tether controller errors to tether commands');
       CsAllocationMat = SIM.parameter('Value',[-1 0 0; 1 0 0; 0 -1 0; 0 0 1],...
           'Unit','','Description','Matrix that allocates control surface controller errors to CS deflections');
    end
    
    methods
        %% constructor
        function obj = threeThrCtlr
            %THREETHRCTLR Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTethers   = SIM.parameter('Description','Number of tethers');
            % altitude tether control
            obj.altiTetherKp   = SIM.parameter('Unit','(m/s)/m','Description','Altitude tether controller Kp');
            obj.altiTetherKi   = SIM.parameter('Unit','(m/s)/(m*s)','Description','Altitude tether controller Ki');
            obj.altiTetherKd   = SIM.parameter('Unit','(m/s)/(m/s)','Description','Altitude tether controller Kd');
            obj.altiTetherTau  = SIM.parameter('Unit','s','Description','Altitude tether controller time constant');
            obj.altiErrorSat   = SIM.parameter('Unit','m','Description','Altitude tether controller error saturation');
            % pitch tether control
            obj.pitchTetherKp  = SIM.parameter('Unit','(m/s)/rad','Description','Pitch tether controller Kp');
            obj.pitchTetherKi  = SIM.parameter('Unit','(m/s)/(rad*s)','Description','Pitch tether controller Ki');
            obj.pitchTetherKd  = SIM.parameter('Unit','(m/s)/(rad/s)','Description','Pitch tether controller Kd');
            obj.pitchTetherTau = SIM.parameter('Unit','s','Description','Pitch tether controller time constant');
            % roll tether control
            obj.rollTetherKp   = SIM.parameter('Unit','(m/s)/rad','Description','Roll tether controller Kp');
            obj.rollTetherKi   = SIM.parameter('Unit','(m/s)/(rad*s)','Description','Roll tether controller Kp');
            obj.rollTetherKd   = SIM.parameter('Unit','(m/s)/(rad/s)','Description','Roll tether controller Kp');
            obj.rollTetherTau  = SIM.parameter('Unit','s','Description','Roll tether controller time constant');
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
        
        % set tether altitude gains
        function setAltiTetherKp(obj,val,units)
            obj.altiTetherKp.setValue(val,units);
        end
        
        function setAltiTetherKi(obj,val,units)
            obj.altiTetherKi.setValue(val,units);
        end

        function setAltiTetherKd(obj,val,units)
            obj.altiTetherKd.setValue(val,units);
        end

        function setAltiTetherTau(obj,val,units)
            obj.altiTetherTau.setValue(val,units);
        end

        function setAltiErrorSat(obj,val,units)
            obj.altiErrorSat.setValue(val,units);
        end
        
        % set tether pitch gains
        function setPitchTetherKp(obj,val,units)
            obj.pitchTetherKp.setValue(val,units);
        end
        
        function setPitchTetherKi(obj,val,units)
            obj.pitchTetherKi.setValue(val,units);
        end
        
        function setPitchTetherKd(obj,val,units)
            obj.pitchTetherKd.setValue(val,units);
        end
        
        function setPitchTetherTau(obj,val,units)
            obj.pitchTetherTau.setValue(val,units);
        end
        
        % set tether roll gains
        function setRollTetherKp(obj,val,units)
            obj.rollTetherKp.setValue(val,units);
        end
        
        function setRollTetherKi(obj,val,units)
            obj.rollTetherKi.setValue(val,units);
        end
        
        function setRollTetherKd(obj,val,units)
            obj.rollTetherKd.setValue(val,units);
        end
        
        function setRollTetherTau(obj,val,units)
            obj.rollTetherTau.setValue(val,units);
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
        function scaleThreeThrCtlr(obj)
            
            LS = obj.lengthScale.Value;
            
            % scale altitude tether control gains
            obj.setAltiTetherKp(obj.altiTetherKp.Value*(1/LS^0.5),'(m/s)/m');
            obj.setAltiTetherKi(obj.altiTetherKi.Value*(1/LS),'(m/s)/(m*s)');
            obj.setAltiTetherTau(obj.altiTetherTau.Value*LS^0.5,'s');
            obj.setAltiErrorSat(obj.altiErrorSat.Value*LS,'m');
            
            % scale pitch tether control gains
            obj.setPitchTetherKp(obj.pitchTetherKp.Value*LS^0.5,'(m/s)/rad');
            obj.setPitchTetherKd(obj.pitchTetherKd.Value*LS,'(m/s)/(rad/s)');
            obj.setPitchTetherTau(obj.pitchTetherTau.Value*LS^0.5,'s');
            
            % scale roll tether control gains
            obj.setRollTetherKp(obj.rollTetherKp.Value*LS^0.5,'(m/s)/rad');
            obj.setRollTetherKd(obj.rollTetherKd.Value*LS,'(m/s)/(rad/s)');
            obj.setRollTetherTau(obj.rollTetherTau.Value*LS^0.5,'s');
            
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


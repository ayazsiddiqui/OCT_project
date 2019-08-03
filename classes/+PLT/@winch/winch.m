classdef winch
    %WINCH Summary of this class goes here
    
    properties (SetAccess = private)
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
            obj.numTethers  = SIM.parameter('Description','Number of tethers');
            obj.wnchMaxTugSpeed = SIM.parameter('Unit','m/s','Description','Max winch pull in (wind up) speed');
            obj.wnchMaxReleaseSpeed = SIM.parameter('Unit','m/s','Description','Max winch push out (wind out) speed');
            obj.wnchTimeConstant = SIM.parameter('Unit','s','Description','Winch time constant');
            obj.initThrLength = SIM.parameter('Unit','m','Description','Initial tether length');
            
        end
        
        %% setters
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
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
        % initial tether length
        function val = recommendInitTetherLength(obj,vehicle,gndStn,tethers,environment)
            
            % calculate total external forces except tethers
            F_grav = vehicle.mass.Value*environment.gravAccel.Value*[0;0;-1];
            F_buoy =  environment.fluidDensity.Value*vehicle.volume.Value*...
                environment.gravAccel.Value*[0;0;1];
            
            % calculate lift forces for wing and HS, ignore VS
            Vrel = environment.inertialFlowVel.Value - vehicle.init_inertialCmVel.Value;
            q = 0.5*environment.fluidDensity.Value*(norm(Vrel))^2;
            Sref = vehicle.fluidRefArea.Value;
            F_aero = [0;0;0];
            
            for ii = 1:3
                CL(ii) = interp1(vehicle.fluidCoeffData(ii).alpha,...
                    vehicle.fluidCoeffData(ii).CL,...
                    (180/pi)*vehicle.init_euler.Value(2));
                CD(ii) = interp1(vehicle.fluidCoeffData(ii).alpha,...
                    vehicle.fluidCoeffData(ii).CD,...
                    (180/pi)*vehicle.init_euler.Value(2));
                F_aero = F_aero + q*Sref*[CD(ii);0;CL(ii)];
            end
            
            sum_F = norm(F_grav + F_buoy + F_aero);
            
            [oCb,~] = rotation_sequence(vehicle.init_euler.Value);
            [oCp,~] = rotation_sequence([0 0 gndStn.init_euler.Value]);
            
            % determine initial tether lenghts
            switch obj.numTethers.Value
                case 1
                    L = norm( vehicle.init_inertialCmPos.Value + ...
                        (oCb*vehicle.thrAttchPts.Value) - ...
                        (oCp*gndStn.thrAttchPts.Value) );
                    
                    delta_L = sum_F/(L*tethers.thrYoungs.Value*...
                        (pi/4)*tethers.thrDiameter.Value^2);
                    
                    val = (L - delta_L);
                    
                case 3
                    L1 = norm( vehicle.init_inertialCmPos.Value + ...
                        (oCb*vehicle.thrAttchPts.Value(:,1)) - ...
                        (oCp*gndStn.thrAttchPts.Value(:,1)) );
                    
                    delta_L1 = (sum_F/4)/(L1*tethers.thrYoungs.Value(1)*...
                        (pi/4)*tethers.thrDiameter.Value(1)^2);
                    
                    % winch 2
                    L2 = norm( vehicle.init_inertialCmPos.Value + ...
                        (oCb*vehicle.thrAttchPts.Value(:,2)) - ...
                        (oCp*gndStn.thrAttchPts.Value(:,2)) );
                    
                    delta_L2 = (sum_F/2)/(L2*tethers.thrYoungs.Value(2)*...
                        (pi/4)*tethers.thrDiameter.Value(2)^2);
                    
                    % winch 3
                    L3 = norm( vehicle.init_inertialCmPos.Value + ...
                        (oCb*vehicle.thrAttchPts.Value(:,3)) - ...
                        (oCp*gndStn.thrAttchPts.Value(:,3)) );
                    
                    delta_L3 = (sum_F/4)/(L3*tethers.thrYoungs.Value(3)*...
                        (pi/4)*tethers.thrDiameter.Value(3)^2);
                    
                    ini_L = [(L1-delta_L1),(L2-delta_L2),(L3-delta_L3)];
                    
                    val = ini_L;
                    
                otherwise
                    error(['Method not progerammed for %d winches.',obj.numTethers])
            end
            
        end
        
    end
end


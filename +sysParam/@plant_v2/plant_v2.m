classdef plant_v2
    %PLANT Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties
    properties (Access = public)
        ScaleFactor
        numTethers
        numTurbines
        vehicle
        turbines
        tethers
        winches
        gndStation
        aeroDataFileName
    end
    
    properties (Dependent)
        aeroCoeffData
        aeroDesignData
        vehicleTetherAttchPts
        gndStationTetherAttchPts
    end
    
    methods
        %% Constructor
        function thisPlant = plant_v2(arg1,arg2)
            %PLANT Construct an instance of this class
            %   Detailed explanation goes here
            
            thisPlant.ScaleFactor = 1;
            % vehicle constants
            thisPlant.vehicle.mass.value = [];
            thisPlant.vehicle.added_mass.value = [];
            thisPlant.vehicle.MI.value = [];
            thisPlant.vehicle.volume.value = [];
            thisPlant.vehicle.Rcb_cm.value = [];
            thisPlant.vehicle.Rcm_wingLE.value = [];
            % vehicle initial conditions
            thisPlant.vehicle.ini_Rcm_o.value = [];
            thisPlant.vehicle.ini_O_Vcm_o.value = [];
            thisPlant.vehicle.ini_euler.value = [];
            thisPlant.vehicle.ini_OwB.value = [];
            % ground station constants
            thisPlant.gndStation.rotationSwitch.value = [];
            thisPlant.gndStation.Izz.value = [];
            thisPlant.gndStation.dampCoeff.value = [];
            % groundstation initial conditions
            thisPlant.gndStation.ini_platform_ang.value = [];
            thisPlant.gndStation.ini_platform_vel.value = [];
            
            if nargin == 1
                thisPlant.numTethers = arg1;
                thisPlant.numTurbines = 2;
                
                fprintf(['Arguments for number of tethers given, setting default value for number of turbines',...
                    'Number of tethers: %d\n Number of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);

                % winches
                thisPlant.winches = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.winches(ii).maxSpeed = [];
                    thisPlant.winches(ii).timeConstant = [];
                    thisPlant.winches(ii).initTetherLength = [];
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.tethers(ii).numNodes = [];
                    thisPlant.tethers(ii).diameter = [];
                    thisPlant.tethers(ii).youngsModulus = [];
                    thisPlant.tethers(ii).dampingRatio = [];
                    thisPlant.tethers(ii).dragCoeff = [];
                    thisPlant.tethers(ii).density = [];
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % turbines
                thisPlant.turbines = struct.empty(0,2);
                for ii = 1:2
                    thisPlant.turbines(ii).Rturb_cm.value = [];
                    thisPlant.turbines(ii).diameter.value = [];
                    thisPlant.turbines(ii).powerCoeff.value = [];
                    thisPlant.turbines(ii).dragCoeff.value = [];
                end
                thisPlant.turbines = reshape(thisPlant.turbines,1,[]);
                
            elseif nargin == 2
                thisPlant.numTethers = arg1;
                thisPlant.numTurbines = arg2;
                
                fprintf(['Arguments for number of tethers and turbines given.\n',...
                    'Number of tethers: %d \nNumber of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);
                
                % winches
                thisPlant.winches = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.winches(ii).maxSpeed = [];
                    thisPlant.winches(ii).timeConstant = [];
                    thisPlant.winches(ii).initTetherLength = [];
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.tethers(ii).numNodes = [];
                    thisPlant.tethers(ii).diameter = [];
                    thisPlant.tethers(ii).youngsModulus = [];
                    thisPlant.tethers(ii).dampingRatio = [];
                    thisPlant.tethers(ii).dragCoeff = [];
                    thisPlant.tethers(ii).density = [];
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % turbines
                thisPlant.turbines = struct.empty(0,arg2);
                for ii = 1:arg2
                    thisPlant.turbines(ii).Rturb_cm.value = [];
                    thisPlant.turbines(ii).diameter.value = [];
                    thisPlant.turbines(ii).powerCoeff.value = [];
                    thisPlant.turbines(ii).dragCoeff.value = [];
                end
                thisPlant.turbines = reshape(thisPlant.turbines,1,[]);
                
            else
                thisPlant.numTethers = 1;
                thisPlant.numTurbines = 2;
                
                fprintf(['Number of tethers or turbines not specified, constructing default plant.\n'...
                    'Number of tethers: %d\n Number of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);
                
                % winches
                thisPlant.winches = struct.empty(0,1);
                for ii = 1:1
                    thisPlant.winches(ii).maxSpeed = [];
                    thisPlant.winches(ii).timeConstant = [];
                    thisPlant.winches(ii).initTetherLength = [];
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct.empty(0,1);
                for ii = 1:1
                    thisPlant.tethers(ii).numNodes = [];
                    thisPlant.tethers(ii).diameter = [];
                    thisPlant.tethers(ii).youngsModulus = [];
                    thisPlant.tethers(ii).dampingRatio = [];
                    thisPlant.tethers(ii).dragCoeff = [];
                    thisPlant.tethers(ii).density = [];
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);

                % turbines
                thisPlant.turbines = struct.empty(0,2);
                for ii = 1:2
                    thisPlant.turbines(ii).Rturb_cm.value = [];
                    thisPlant.turbines(ii).diameter.value = [];
                    thisPlant.turbines(ii).powerCoeff.value = [];
                    thisPlant.turbines(ii).dragCoeff.value = [];
                end
                thisPlant.turbines = reshape(thisPlant.turbines,1,[]);
            end
            
        end
        
        %% Get methods
        
        function val = get.aeroCoeffData(obj)
            val = load(obj.aeroDataFileName,'aeroStruct');
            val = val.aeroStruct;
            for ii = 1:4
                val(ii).aeroCentPosVec = -obj.vehicle.Rcm_wingLE.value...
                    + val(ii).aeroCentPosVec;
            end
        end
        
        function val = get.aeroDesignData(obj)
            val = load(obj.aeroDataFileName,'dsgnData');
            val = val.dsgnData;
        end
        
        function val = get.vehicleTetherAttchPts(obj)
            
            val = struct.empty(0,obj.numTethers);
            if obj.numTethers == 1
                val(1).value = [0;0;0];
                
            elseif obj.numTethers == 3
                val = struct.empty(0,3);
                % port tether
                val(1).value = -obj.vehicle.Rcm_wingLE.value + ...
                    [(tand(obj.aeroDesignData.wing_sweep)*obj.aeroDesignData.wing_span/2) + ...
                    obj.aeroDesignData.wing_chord*obj.aeroDesignData.wing_TR/4;...
                    -obj.aeroDesignData.wing_span/2; ...
                    tand(obj.aeroDesignData.wing_dihedral)*obj.aeroDesignData.wing_span/2];
                % aft tether
                val(2).value = -obj.vehicle.Rcm_wingLE.value + ...
                    [obj.aeroDesignData.h_stab_LE + (obj.aeroDesignData.h_stab_chord);0 ;0];
                % starboard tether
                val(3).value = -obj.vehicle.Rcm_wingLE.value + ...
                    [(tand(obj.aeroDesignData.wing_sweep)*obj.aeroDesignData.wing_span/2) + ...
                    obj.aeroDesignData.wing_chord*obj.aeroDesignData.wing_TR/4;...
                    obj.aeroDesignData.wing_span/2; ...
                    tand(obj.aeroDesignData.wing_dihedral)*obj.aeroDesignData.wing_span/2];
            else
                error('Method not defined for %d tethers',obj.numTethers);
                
            end
        end
        
        function val = get.gndStationTetherAttchPts(obj)
            
            val = struct.empty(0,obj.numTethers);
            for ii = 1:obj.numTethers
                val(ii).value = obj.vehicleTetherAttchPts(ii).value;
                val(ii).value(3) = 0;
            end
        end
        
        %% set methods
        
        
        %% call methods
        % design tether diameter based on max flow speed and aloowable
        % elongation
        function obj = designTetherDiameter(obj,env,maxAppFlowMultiplier,maxPercentageElongation)
            % calculate total external forces except tethers
            F_grav = obj.vehicle.mass.value*env.gravAccel.value*[0;0;-1];
            F_buoy =  env.flowDensity.value*obj.vehicle.volume.value*...
                env.gravAccel.value*[0;0;1];
            
            % calculate lift forces for wing and HS, ignore VS
            q_max = 0.5*env.flowDensity.value*(maxAppFlowMultiplier*norm(env.iniertialFlowVel.value))^2;
            Sref = obj.aeroCoeffData(1).refArea;
            F_aero = [0;0;0];
            for ii = 1:3
                CLm(ii) = max(obj.aeroCoeffData(ii).CL);
                F_aero = F_aero + q_max*Sref*[0;0;CLm(ii)];
            end
            
            sum_F = norm(F_grav + F_buoy + F_aero);
            
            switch obj.numTethers
                case 1
                    obj.tethers(1).diameter = sqrt((4*sum_F)/...
                        (pi*maxPercentageElongation*obj.tethers(1).youngsModulus));
                case 3
                    obj.tethers(1).diameter = sqrt((4*sum_F/4)/...
                        (pi*maxPercentageElongation*obj.tethers(1).youngsModulus));
                    obj.tethers(2).diameter = sqrt((4*sum_F/2)/...
                        (pi*maxPercentageElongation*obj.tethers(2).youngsModulus));
                    obj.tethers(3).diameter = sqrt((4*sum_F/4)/...
                        (pi*maxPercentageElongation*obj.tethers(3).youngsModulus));
                otherwise
                    error(['What are you trying to achieve by running this system with %d tether?! '...
                        'I didn''t account for that!\n',obj.numTethers])
            end
        end
        
        % calculate initial tether unstretched length
        function obj = setTetherInitLength(obj,env)
            
            % calculate total external forces except tethers
            F_grav = obj.vehicle.mass.value*env.gravAccel.value*[0;0;-1];
            F_buoy =  env.flowDensity.value*obj.vehicle.volume.value*...
                env.gravAccel.value*[0;0;1];
            
            % calculate lift forces for wing and HS, ignore VS
            Vrel = env.iniertialFlowVel.value - obj.vehicle.ini_O_Vcm_o.value;
            q = 0.5*env.flowDensity.value*(norm(Vrel))^2;
            Sref = obj.aeroCoeffData(1).refArea;
            F_aero = [0;0;0];
            
            for ii = 1:3
                CL(ii) = interp1(obj.aeroCoeffData(ii).alpha,...
                    obj.aeroCoeffData(ii).CL,...
                    (180/pi)*obj.vehicle.ini_euler.value(2));
                CD(ii) = interp1(obj.aeroCoeffData(ii).alpha,...
                    obj.aeroCoeffData(ii).CD,...
                    (180/pi)*obj.vehicle.ini_euler.value(2));
                F_aero = F_aero + q*Sref*[CD(ii);0;CL(ii)];
            end
            
            sum_F = norm(F_grav + F_buoy + F_aero);
            
            [oCb,~] = rotation_sequence(obj.vehicle.ini_euler.value);
            [oCp,~] = rotation_sequence([0 0 obj.gndStation.ini_platform_ang.value]);
            
            switch obj.numTethers
                case 1
                    L = norm( obj.vehicle.ini_Rcm_o.value + ...
                        (oCb*obj.vehicleTetherAttchPts(1).value) - ...
                        (oCp*obj.gndStationTetherAttchPts(1).value) );
                    
                    delta_L = sum_F/(L*obj.tethers(1).youngsModulus*...
                        (pi/4)*obj.tethers(1).diameter^2);
                    
                    obj.winches(1).initTetherLength = (L - delta_L);
                    
                case 3
                    L1 = norm( obj.vehicle.ini_Rcm_o.value + ...
                        (oCb*obj.vehicleTetherAttchPts(1).value) - ...
                        (oCp*obj.gndStationTetherAttchPts(1).value) );
                    
                    delta_L1 = (sum_F/4)/(L1*obj.tethers(1).youngsModulus*...
                        (pi/4)*obj.tethers(1).diameter^2);
                    
                    obj.winches(1).initTetherLength = (L1 - delta_L1);
                    
                    % winch 2
                    L2 = norm( obj.vehicle.ini_Rcm_o.value + ...
                        (oCb*obj.vehicleTetherAttchPts(2).value) - ...
                        (oCp*obj.gndStationTetherAttchPts(2).value) );
                    
                    delta_L2 = (sum_F/2)/(L2*obj.tethers(2).youngsModulus*...
                        (pi/4)*obj.tethers(2).diameter^2);
                    
                    obj.winches(2).initTetherLength = (L2 - delta_L2);
                    
                    % winch 3
                    L3 = norm( obj.vehicle.ini_Rcm_o.value + ...
                        (oCb*obj.vehicleTetherAttchPts(3).value) - ...
                        (oCp*obj.gndStationTetherAttchPts(3).value) );
                    
                    delta_L3 = (sum_F/2)/(L3*obj.tethers(3).youngsModulus*...
                        (pi/4)*obj.tethers(3).diameter^2);
                    
                    obj.winches(3).initTetherLength = (L3 - delta_L3);
                    
                otherwise
                    error(['Method not progerammed for %d winches.',obj.numTethers])
            end
            
        end
        
        % calc added mass
        function obj = calcAddedMass(obj,env)
            density = env.flowDensity.value;
            span = obj.aeroDesignData.wing_span;
            chord = obj.aeroDesignData.wing_chord;
            HS_span = obj.aeroDesignData.h_stab_span;
            HS_chord = obj.aeroDesignData.h_stab_chord;
            VS_span = obj.aeroDesignData.v_stab_span;
            VS_chord = obj.aeroDesignData.v_stab_chord;
            
            m_added_x = pi*density*(span*(0.15*chord/2)^2 + ...
                HS_span*(0.15*HS_chord/2)^2 + VS_span*(0.15*VS_chord/2)^2);
            m_added_y = pi*density*(1.98*span*(chord/2)^2 + ...
                1.98*HS_span*(HS_chord/2)^2 + VS_span*(VS_chord/2)^2);
            m_added_z = pi*density*(span*(chord/2)^2 + ...
                HS_span*(HS_chord/2)^2 + 1.98*VS_span*(VS_chord/2)^2);
            
            obj.vehicle.added_mass.value = [m_added_x 0 0;0 m_added_y 0; 0 0 m_added_z];
        end
        
        
    end
end


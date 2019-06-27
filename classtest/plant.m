classdef plant
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
    end
    
    properties (Dependent)
        added_mass
    end
    
    methods
        %% Constructor
        function thisPlant = plant(arg1,arg2)
            %PLANT Construct an instance of this class
            %   Detailed explanation goes here
            thisPlant.ScaleFactor = 1;
            thisPlant.vehicle.mass = param([],'kg','Mass of the lifting body');
            thisPlant.vehicle.MI = param(eye(3),'kg-m^2','Moment of inertia about CM in body frame');
            thisPlant.vehicle.volume = param([],'m^3','Volume of lifting body');
            thisPlant.vehicle.Rcb_cm = param([],'m','Vector going from CM to center of buoyancy');
            thisPlant.vehicle.Rcm_wingLE = param([],'m','Vector going from wing LE to CM');
            thisPlant.gndStation.Izz = param([],'kg-m^2','Moment of inertia about the ground frame Z axis');
            thisPlant.gndStation.dampCoeff = param([],'N-m-s','Damping coeffient of platform');
            
            if nargin == 1
                thisPlant.numTethers = arg1;
                thisPlant.numTurbines = 2;
                
                fprintf(['Arguments for number of tethers given, setting default value for number of turbines',...
                    'Number of tethers: %d\n Number of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);

                % vehicle tether attachment points
                thisPlant.vehicle.tetherAttchPts = param.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.vehicle.tetherAttchPts(ii) = param([],'m',...
                        'Vector going from vehicle CM to tether attachment point');
                end
                thisPlant.vehicle.tetherAttchPts = reshape(thisPlant.vehicle.tetherAttchPts,1,[]);
                
                % winches
                thisPlant.winches = struct([]);
                for ii = 1:arg1
                    thisPlant.winches(ii).maxSpeed = param([],'m/s','Maximum spooling speed');
                    thisPlant.winches(ii).timeConstant = param([],'s','Winch time constant');
                    thisPlant.winches(ii).initTetherLength = param([],'m',...
                        'Initial unstretched tether lengths');
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct([]);
                for ii = 1:arg1
                    thisPlant.tethers(ii).numNodes = param([],'','Number of nodes');
                    thisPlant.tethers(ii).diameter = param([],'m','Tether diameter');
                    thisPlant.tethers(ii).youngsModulus = param([],'N/m^2',...
                        'Tether Youngs modulus');
                    thisPlant.tethers(ii).dampingRation = param([],'',...
                        'Tether damping ratio');
                    thisPlant.tethers(ii).dragCoeff = param([],'',...
                        'Tether drag Coefficient');
                    thisPlant.tethers(ii).density = param([],'kg/m^3',...
                        'Tether density');
                    thisPlant.tethers(ii).vehicleMass = thisPlant.vehicle.mass;
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % ground station attachment points
                thisPlant.gndStation.tetherAttchPts = param.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.gndStation.tetherAttchPts(ii) = param([],'m',...
                        'Vector going from vehicle CM to tether attachment point');
                end
                thisPlant.gndStation.tetherAttchPts = reshape(thisPlant.gndStation.tetherAttchPts,1,[]);

                % turbines
                thisPlant.turbines = struct([]);
                for ii = 1:2
                    thisPlant.turbines(ii).Rturb_cm = param([],'m',...
                        'Vector going from vehicle CM to turbine');
                    thisPlant.turbines(ii).diameter = param([],'m',...
                        'Turbine diameter');
                    thisPlant.turbines(ii).powerCoeff = param([],'',...
                        'Turbine power coefficient');
                    thisPlant.turbines(ii).dragCoeff = param([],'',...
                        'Turbine drag coefficient');
                end
                thisPlant.turbines = reshape(thisPlant.turbines,1,[]);
                
            elseif nargin == 2
                thisPlant.numTethers = arg1;
                thisPlant.numTurbines = arg2;
                
                fprintf(['Arguments for number of tethers and turbines given.\n',...
                    'Number of tethers: %d \nNumber of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);
                
                thisPlant.vehicle.tetherAttchPts = param.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.vehicle.tetherAttchPts(ii) = param([],'m',...
                        'Vector going from vehicle CM to tether attachment point');
                end
                thisPlant.vehicle.tetherAttchPts = reshape(thisPlant.vehicle.tetherAttchPts,1,[]);
                
                % winches
                thisPlant.winches = struct([]);
                for ii = 1:arg1
                    thisPlant.winches(ii).maxSpeed = param([],'m/s','Maximum spooling speed');
                    thisPlant.winches(ii).timeConstant = param([],'s','Winch time constant');
                    thisPlant.winches(ii).initTetherLength = param([],'m',...
                        'Initial unstretched tether lengths');
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct([]);
                for ii = 1:arg1
                    thisPlant.tethers(ii).numNodes = param([],'','Number of nodes');
                    thisPlant.tethers(ii).diameter = param([],'m','Tether diameter');
                    thisPlant.tethers(ii).youngsModulus = param([],'N/m^2',...
                        'Tether Youngs modulus');
                    thisPlant.tethers(ii).dampingRation = param([],'',...
                        'Tether damping ratio');
                    thisPlant.tethers(ii).dragCoeff = param([],'',...
                        'Tether drag Coefficient');
                    thisPlant.tethers(ii).density = param([],'kg/m^3',...
                        'Tether density');
                    thisPlant.tethers(ii).vehicleMass = thisPlant.vehicle.mass;
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % ground station attachment points
                thisPlant.gndStation.tetherAttchPts = param.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.gndStation.tetherAttchPts(ii) = param([],'m',...
                        'Vector going from vehicle CM to tether attachment point');
                end
                thisPlant.gndStation.tetherAttchPts = reshape(thisPlant.gndStation.tetherAttchPts,1,[]);

                %turbines
                thisPlant.turbines = struct([]);
                for ii = 1:arg2
                    thisPlant.turbines(ii).Rturb_cm = param([],'m',...
                        'Vector going from vehicle CM to turbine');
                    thisPlant.turbines(ii).diameter = param([],'m',...
                        'Turbine diameter');
                    thisPlant.turbines(ii).powerCoeff = param([],'',...
                        'Turbine power coefficient');
                    thisPlant.turbines(ii).dragCoeff = param([],'',...
                        'Turbine drag coefficient');
                end
                thisPlant.turbines = reshape(thisPlant.turbines,1,[]);

            else
                thisPlant.numTethers = 1;
                thisPlant.numTurbines = 2;
                
                fprintf(['Number of tethers or turbines not specified, constructing default plant.\n'...
                    'Number of tethers: %d\n Number of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);
                
                thisPlant.vehicle.tetherAttchPts = param.empty(0,1);
                for ii = 1:1
                    thisPlant.vehicle.tetherAttchPts(ii) = param([],'m',...
                        'Vector going from vehicle CM to tether attachment point');
                end
                thisPlant.vehicle.tetherAttchPts = reshape(thisPlant.vehicle.tetherAttchPts,1,[]);
                
                % winches
                thisPlant.winches = struct([]);
                for ii = 1:1
                    thisPlant.winches(ii).maxSpeed = param([],'m/s','Maximum spooling speed');
                    thisPlant.winches(ii).timeConstant = param([],'s','Winch time constant');
                    thisPlant.winches(ii).initTetherLength = param([],'m',...
                        'Initial unstretched tether lengths');
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct([]);
                for ii = 1:1
                    thisPlant.tethers(ii).numNodes = param([],'','Number of nodes');
                    thisPlant.tethers(ii).diameter = param([],'m','Tether diameter');
                    thisPlant.tethers(ii).youngsModulus = param([],'N/m^2',...
                        'Tether Youngs modulus');
                    thisPlant.tethers(ii).dampingRation = param([],'',...
                        'Tether damping ratio');
                    thisPlant.tethers(ii).dragCoeff = param([],'',...
                        'Tether drag Coefficient');
                    thisPlant.tethers(ii).density = param([],'kg/m^3',...
                        'Tether density');
                    thisPlant.tethers(ii).vehicleMass = thisPlant.vehicle.mass;
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % ground station attachment points
                thisPlant.gndStation.tetherAttchPts = param.empty(0,1);
                for ii = 1:1
                    thisPlant.gndStation.tetherAttchPts(ii) = param([],'m',...
                        'Vector going from vehicle CM to tether attachment point');
                end
                thisPlant.gndStation.tetherAttchPts = reshape(thisPlant.gndStation.tetherAttchPts,1,[]);

                %turbines
                thisPlant.turbines = struct([]);
                for ii = 1:2
                    thisPlant.turbines(ii).Rturb_cm = param([],'m',...
                        'Vector going from vehicle CM to turbine');
                    thisPlant.turbines(ii).diameter = param([],'m',...
                        'Turbine diameter');
                    thisPlant.turbines(ii).powerCoeff = param([],'',...
                        'Turbine power coefficient');
                    thisPlant.turbines(ii).dragCoeff = param([],'',...
                        'Turbine drag coefficient');
                end
                thisPlant.turbines = reshape(thisPlant.turbines,1,[]);

            end
            
        end
        
        %% Get methods
        
        %% set methods
        
    end
end


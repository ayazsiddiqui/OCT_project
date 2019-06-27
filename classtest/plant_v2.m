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
    end
    
    properties (Dependent)
        added_mass
    end
    
    methods
        %% Constructor
        function thisPlant = plant_v2(arg1,arg2)
            %PLANT Construct an instance of this class
            %   Detailed explanation goes here
            thisPlant.ScaleFactor = 1;
            thisPlant.vehicle.mass.value = [];
            thisPlant.vehicle.MI.value = [];
            thisPlant.vehicle.volume.value = [];
            thisPlant.vehicle.Rcb_cm.value = [];
            thisPlant.vehicle.Rcm_wingLE.value = [];
            thisPlant.gndStation.Izz.value = [];
            thisPlant.gndStation.dampCoeff.value = [];
            
            if nargin == 1
                thisPlant.numTethers = arg1;
                thisPlant.numTurbines = 2;
                
                fprintf(['Arguments for number of tethers given, setting default value for number of turbines',...
                    'Number of tethers: %d\n Number of turbines: %d\n'],...
                    thisPlant.numTethers,thisPlant.numTurbines);

                % vehicle tether attachment points
                thisPlant.vehicle.tetherAttchPts = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.vehicle.tetherAttchPts(ii).value = [];
                end
                thisPlant.vehicle.tetherAttchPts = reshape(thisPlant.vehicle.tetherAttchPts,1,[]);
                
                % winches
                thisPlant.winches = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.winches(ii).maxSpeed.value = [];
                    thisPlant.winches(ii).timeConstant.value = [];
                    thisPlant.winches(ii).initTetherLength.value = [];
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.tethers(ii).numNodes.value = [];
                    thisPlant.tethers(ii).diameter.value = [];
                    thisPlant.tethers(ii).youngsModulus.value = [];
                    thisPlant.tethers(ii).dampingRation.value = [];
                    thisPlant.tethers(ii).dragCoeff.value = [];
                    thisPlant.tethers(ii).density.value = [];
                    thisPlant.tethers(ii).vehicleMass.value = thisPlant.vehicle.mass;
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % ground station attachment points
                thisPlant.gndStation.tetherAttchPts = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.gndStation.tetherAttchPts(ii).value = [];
                end
                thisPlant.gndStation.tetherAttchPts = reshape(thisPlant.gndStation.tetherAttchPts,1,[]);

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
                
                % vehicle tether attachment points
                thisPlant.vehicle.tetherAttchPts = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.vehicle.tetherAttchPts(ii).value = [];
                end
                thisPlant.vehicle.tetherAttchPts = reshape(thisPlant.vehicle.tetherAttchPts,1,[]);
                
                % winches
                thisPlant.winches = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.winches(ii).maxSpeed.value = [];
                    thisPlant.winches(ii).timeConstant.value = [];
                    thisPlant.winches(ii).initTetherLength.value = [];
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.tethers(ii).numNodes.value = [];
                    thisPlant.tethers(ii).diameter.value = [];
                    thisPlant.tethers(ii).youngsModulus.value = [];
                    thisPlant.tethers(ii).dampingRation.value = [];
                    thisPlant.tethers(ii).dragCoeff.value = [];
                    thisPlant.tethers(ii).density.value = [];
                    thisPlant.tethers(ii).vehicleMass.value = thisPlant.vehicle.mass;
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % ground station attachment points
                thisPlant.gndStation.tetherAttchPts = struct.empty(0,arg1);
                for ii = 1:arg1
                    thisPlant.gndStation.tetherAttchPts(ii).value = [];
                end
                thisPlant.gndStation.tetherAttchPts = reshape(thisPlant.gndStation.tetherAttchPts,1,[]);

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
                
                % vehicle tether attachment points
                thisPlant.vehicle.tetherAttchPts = struct.empty(0,1);
                for ii = 1:1
                    thisPlant.vehicle.tetherAttchPts(ii).value = [];
                end
                thisPlant.vehicle.tetherAttchPts = reshape(thisPlant.vehicle.tetherAttchPts,1,[]);
                
                % winches
                thisPlant.winches = struct.empty(0,1);
                for ii = 1:1
                    thisPlant.winches(ii).maxSpeed.value = [];
                    thisPlant.winches(ii).timeConstant.value = [];
                    thisPlant.winches(ii).initTetherLength.value = [];
                end
                thisPlant.winches = reshape(thisPlant.winches,1,[]);
                
                % tethers
                thisPlant.tethers = struct.empty(0,1);
                for ii = 1:1
                    thisPlant.tethers(ii).numNodes.value = [];
                    thisPlant.tethers(ii).diameter.value = [];
                    thisPlant.tethers(ii).youngsModulus.value = [];
                    thisPlant.tethers(ii).dampingRation.value = [];
                    thisPlant.tethers(ii).dragCoeff.value = [];
                    thisPlant.tethers(ii).density.value = [];
                    thisPlant.tethers(ii).vehicleMass.value = thisPlant.vehicle.mass;
                end
                thisPlant.tethers = reshape(thisPlant.tethers,1,[]);
                
                % ground station attachment points
                thisPlant.gndStation.tetherAttchPts = struct.empty(0,1);
                for ii = 1:1
                    thisPlant.gndStation.tetherAttchPts(ii).value = [];
                end
                thisPlant.gndStation.tetherAttchPts = reshape(thisPlant.gndStation.tetherAttchPts,1,[]);

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
        function val = get.added_mass(obj)
            val = param([1.8017e+04 0 0;0 1.5825e+06 0;0 0 8.0374e+05],...
                'kg','Added mass on vehicle');
        end
        
        %% set methods
        
    end
end


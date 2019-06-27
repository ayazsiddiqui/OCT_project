classdef plant
    %PLANT Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties
    properties (Access = public)
        ScaleFactor
        vehicle
        gndStation
        numTethers
        numTurbines
    end
    
    properties (Dependant)
        vehicleTetherAttchPts
        turbines
        tethers
        winches
        gndStationTetherAttchPts
        added_mass
    end
    
    methods
        %% Constructor
        function thisPlant = plant()
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
            thisPlant.numTethers = 1;
            thisPlant.numTurbines = 2;
        end
        
        %% Get methods
        %%%% intialize vehicle tether attachment points %%%%%
        function val = get.vehicleTetherAttchPts(somePlant)
            for ii = 1:somePlant.numTethers
                val(ii) = param([],'m',...
                    'Vector going from vehicle CM to tether attachment point');
            end
            val = reshape(val,1,[]);
        end
        
        %%%%% initialize turbines %%%%%
        function val = get.turbines(somePlant)
            val(somePlant.numTurbines) = struct();
            for ii = 1:somePlant.numTurbines
                val(ii).Rturb_cm = param([],'m',...
                    'Vector going from vehicle CM to turbine');
                val(ii).diameter = param([],'m',...
                    'Turbine diameter');
                val(ii).powerCoeff = param([],'',...
                    'Turbine power coefficient');
                val(ii).dragCoeff = param([],'',...
                    'Turbine drag coefficient');
            end
            val = reshape(val,1,[]);
        end
        
        %%%%% initialize tethers %%%%%
        function val = get.tethers(somePlant)
            val(somePlant.numTethers) = struct();
            for ii = 1:somePlant.numTethers
                val(ii).numNodes = param([],'','Number of nodes');
                val(ii).diameter = param([],'m','Tether diameter');
                val(ii).youngsModulus = param([],'N/m^2',...
                    'Tether Youngs modulus');
                val(ii).dampingRation = param([],'',...
                    'Tether damping ratio');
                val(ii).dragCoeff = param([],'',...
                    'Tether drag Coefficient');
                val(ii).density = param([],'kg/m^3',...
                    'Tether density');
                val(ii).vehicleMass = somePlant.vehicle.mass;
            end
            val = reshape(val,1,[]);
        end
        
        %%%%% initialize winches %%%%%
        function val = get.winches(somePlant)
            val(somePlant.numTethers) = struct();
            for ii = 1:somePlant.numTethers
                val(ii).maxSpeed = param([],'m/s','Maximum spooling speed');
                val(ii).timeConstant = param([],'s','Winch time constant');
                val(ii).initTetherLength = param([],'m',...
                    'Initial unstretched tether lengths');
            end
            val = reshape(val,1,[]);
        end
        
        %%%%% initialize ground station tether attachment points %%%%%
        function val = get.gndStationTetherAttchPts(somePlant)
            for ii = 1:somePlant.numTethers
                val(ii) = param([],'m',...
                    'Vector going from ground station CM to tether attachment point');
            end
            val = reshape(val,1,[]);
        end
        
        %%%%% calculate added mass
        function val = get.added_mass(somePlant)
        val = param([1.8017e+04 0 0;0 1.5825e+06 0;0 0 8.0374e+05],...
            'm','Added mass matrix');
            
        end
        
        %% set methods
        %%%%% vehicle tether attachment points %%%%%
        function somePlant = set.vehicleTetherAttchPts(somePlant,value)
            for ii = 1:1:somePlant.numTurbines
            somePlant.vehicleTetherAttchPts(ii).value = value;
            end
        end
        
    end
end


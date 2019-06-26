classdef vehicle 
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lengthScaleFactor = 1   % master length scale factor
        Rcb_cm                  % vec from CM to center of buoy                    
        Rcm_wingLE              % vec from wing LE to CM
        volume                  % body volume
        BF                      % factor of buoyancy
        numTurbine              % number of turbines
        numTether               % number of tethers
        
    end
    
    properties (Dependent)
        flu
    end
    
    methods
        function obj = vehicle
            %VEHICLE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Rcb_cm = systemParams.param('Unit','m');
            obj.Rcm_wingLE = systemParams.param('Unit','m');
            obj.volume = systemParams.param('Unit','m^3');
            
        end

        
        function obj = Scale(obj)
            %Scale parameters
            obj.Rcb_cm.Value = obj.Rcb_cm.Value/obj.lengthScaleFactor;
        end
        
        function val = get.flu(obj)
            val = 2*obj.Rcb_cm.Value;
        end
        
        
    end
end


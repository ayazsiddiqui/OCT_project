classdef vehicle < dynamicprops
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lengthScaleFactor   % master length scale factor
        Rcb_cm              % vec from CM to center of buoy                    
        Rcm_wingLE          % vec from wing LE to CM
        volume              % body volume
        
        
    end
    
    methods
%         function obj = vehicle(inputArg1,inputArg2)
%             %VEHICLE Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


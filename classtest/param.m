classdef param
    
    properties (Access = public)
        value
        units
        description
    end
    
    methods
        function thisParam = param(value, units, description)
            if nargin == 3
                thisParam.value = value;
                thisParam.units = units;
                thisParam.description = description;
            elseif nargin == 2
                thisParam.value = value;
                thisParam.units = units;
                thisParam.description = char();
            elseif nargin == 1
                thisParam.value = value;
                thisParam.units = char();
                thisParam.description = char();
            else
                thisParam.value = [];
                thisParam.units = char();
                thisParam.description = char();
            end
            
        end
    end
end


classdef param < Simulink.Parameter
    
    methods
        function obj = param(varargin)
            p = inputParser;
            addParameter(p,'Value',[],@isnumeric)
            addParameter(p,'Min',[],@isnumeric)
            addParameter(p,'Max',[],@isnumeric)
            addParameter(p,'Unit','',@ischar)
            addParameter(p,'Description','',@ischar)
            parse(p,varargin{:})
            
            for ii = 1:length(p.Parameters)
                obj.(p.Parameters{ii}) = p.Results.(p.Parameters{ii});
            end
        end
        
    end
end


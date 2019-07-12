classdef tether
    %TETHER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        numNodes
        thrDiameter
        thrDensity
        thrYoungs
        thrDampingRatio
        thrDragCoeff
        
    end
    
    methods
        %% constructor
        function obj = tether
            %TETHER Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTethers  = SIM.parameter('Description','Number of tethers');
            obj.numNodes = SIM.parameter('Description','Number of nodes on each tether');
            obj.thrDiameter = SIM.parameter('Unit','m','Description','Tether Diameter');
            obj.thrDensity = SIM.parameter('Unit','kg/m^3','Description','Tether density');
            obj.thrYoungs = SIM.parameter('Unit','N/m^2','Description','Tether Young''s modulus');
            obj.thrDampingRatio = SIM.parameter('Description','Tether damping ratio');
            obj.thrDragCoeff = SIM.parameter('Description','Tether drag coefficient');
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
        
        function setNumNodes(obj,val,units)
            if numel(1) ~= 1
                error('Please provide a single value for the number of nodes');
            else
                obj.numNodes.setValue(val,units);
            end
        end
        
        function setThrDiameter(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.thrDiameter.setValue(reshape(val,1,[]),units);
            end
        end
        
        function setThrDensity(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.thrDensity.setValue(reshape(val,1,[]),units);
            end
        end
        
        function setThrYoungs(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.thrYoungs.setValue(reshape(val,1,[]),units);
            end
        end
        
        function setThrDampingRatio(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.thrDampingRatio.setValue(reshape(val,1,[]),units);
            end
        end
        
        function setThrDragCoeff(obj,val,units)
            if numel(val) ~= obj.numTethers.Value
                error('Number of values provided not equal to number of tethers');
            else
                obj.thrDragCoeff.setValue(reshape(val,1,[]),units);
            end
        end
        
        %% other methods
        
        % scale tethers
        function scaleTether(obj)
            LS = obj.lengthScale.Value;
            DS = obj.densityScale.Value;
            
            obj.setThrDiameter(obj.thrDiameter.Value.*(LS*DS^(1/2)),'m');
            obj.setThrDensity(obj.thrDensity.Value.*DS,'kg/m^3');
            obj.setThrYoungs(obj.thrYoungs.Value.*LS,'N/m^2');
        end
        
        
            
    end
end


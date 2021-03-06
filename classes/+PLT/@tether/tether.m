classdef tether
    %TETHER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
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
            obj.numTethers  = SIM.parameter('Description','Number of tethers');
            obj.numNodes = SIM.parameter('Description','Number of nodes on each tether');
            obj.thrDiameter = SIM.parameter('Unit','m','Description','Tether Diameter');
            obj.thrDensity = SIM.parameter('Unit','kg/m^3','Description','Tether density');
            obj.thrYoungs = SIM.parameter('Unit','N/m^2','Description','Tether Young''s modulus');
            obj.thrDampingRatio = SIM.parameter('Description','Tether damping ratio');
            obj.thrDragCoeff = SIM.parameter('Description','Tether drag coefficient');
        end
        
        %% setters
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
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
        % design tether diameter
        function val = recommendTetherDiameter...
                (obj,vehicle,environment,maxAppFlowMultiplier,maxPercentageElongation)
            
            % calculate total external forces except tethers
            F_grav = vehicle.mass.Value*environment.gravAccel.Value*[0;0;-1];
            F_buoy =  environment.fluidDensity.Value*vehicle.volume.Value*...
                environment.gravAccel.Value*[0;0;1];
            
            % calculate lift forces for wing and HS, ignore VS
            q_max = 0.5*environment.fluidDensity.Value...
                *(maxAppFlowMultiplier*norm(environment.inertialFlowVel.Value))^2;
            Sref = vehicle.fluidRefArea.Value;
            F_aero = [0;0;0];
            for ii = 1:3
                CLm(ii) = max(vehicle.fluidCoeffData(ii).CL);
                F_aero = F_aero + q_max*Sref*[0;0;CLm(ii)];
            end
            
            sum_F = norm(F_grav + F_buoy + F_aero);
            
            switch obj.numTethers.Value
                case 1
                    td1 = sqrt((4*sum_F)/...
                        (pi*maxPercentageElongation*obj.thrYoungs.Value));
                    val = td1;
                    
                case 3
                    td1 = sqrt((4*sum_F/4)/...
                        (pi*maxPercentageElongation*obj.thrYoungs.Value(1)));
                    td2 = sqrt((4*sum_F/2)/...
                        (pi*maxPercentageElongation*obj.thrYoungs.Value(2)));
                    td3 = sqrt((4*sum_F/4)/...
                        (pi*maxPercentageElongation*obj.thrYoungs.Value(3)));
                    
                    val = [td1,td2,td3];
                    
                otherwise
                    error(['What are you trying to achieve by running this system with %d tether?! '...
                        'I didn''t account for that!\n',obj.numTethers])
            end
            
        end
        
        
        
    end
end


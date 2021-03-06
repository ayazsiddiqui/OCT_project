classdef vehicle < dynamicprops
    
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        numTethers
        numTurbines
        buoyFactor
        fluidDensity
        volume
        Ixx
        Iyy
        Izz
        Ixy
        Ixz
        Iyz
        Rcb_cm
        Rbridle_cm
        % aero properties
        % data file name
        fluidCoeffsFileName
        % wing
        RwingLE_cm
        wingChord
        wingAR
        wingTR
        wingSweep
        wingDihedral
        wingIncidence
        wingNACA
        wingClMax
        wingClMin
        % H-stab
        RhsLE_wingLE
        hsChord
        hsAR
        hsTR
        hsSweep
        hsDihedral
        hsIncidence
        hsNACA
        hsClMax
        hsClMin
        % V-stab
        Rvs_wingLE
        vsChord
        vsSpan
        vsTR
        vsSweep
        vsNACA
        vsClMax
        vsClMin
        % intial conditions
        init_inertialCmPos
        init_inertialCmVel
        init_euler
        init_angVel
    end
    
    properties
        % fluid dynamic data
        fluidCoeffData
    end
    
    properties (Dependent)
        mass
        MI
        addedMass
        surfaceOutlines
        thrAttchPts
        turbineAttchPts
        fluidMomentArms
        fluidRefArea
    end
    
    methods
        %% Constructor
        function obj = vehicle
            %VEHICLE Construct an instance of this class
            obj.numTethers  = SIM.parameter('Description','Number of tethers');
            obj.numTurbines = SIM.parameter('Description','Number of turbines');
            obj.buoyFactor = SIM.parameter('Description','Buoyancy Factor');
            obj.fluidDensity = SIM.parameter('Unit','kg/m^3','Description','Fluid density');
            % mass, volume and inertia
            obj.volume         = SIM.parameter('Unit','m^3','Description','volume');
            obj.Ixx            = SIM.parameter('Unit','kg*m^2','Description','Ixx');
            obj.Iyy            = SIM.parameter('Unit','kg*m^2','Description','Iyy');
            obj.Izz            = SIM.parameter('Unit','kg*m^2','Description','Izz');
            obj.Ixy            = SIM.parameter('Unit','kg*m^2','Description','Ixy');
            obj.Ixz            = SIM.parameter('Unit','kg*m^2','Description','Ixz');
            obj.Iyz            = SIM.parameter('Unit','kg*m^2','Description','Iyz');
            % some vectors
            obj.Rcb_cm        = SIM.parameter('Unit','m','Description','Vector going from CM to center of buoyancy');
            obj.Rbridle_cm        = SIM.parameter('Unit','m','Description','Vector going from CM to bridle point');
            % fluid coeffs file name
            obj.fluidCoeffsFileName = SIM.parameter('Description','File that contains fluid dynamics coefficient data');
            % 
            obj.fluidCoeffData = [];
            % defining aerodynamic surfaces
            obj.RwingLE_cm    = SIM.parameter('Unit','m','Description','Vector going from CM to wing leading edge');
            obj.wingChord     = SIM.parameter('Unit','m','Description','Wing root chord');
            obj.wingAR        = SIM.parameter('Description','Wing Aspect ratio');
            obj.wingTR        = SIM.parameter('Description','Wing Taper ratio');
            obj.wingSweep     = SIM.parameter('Unit','deg','Description','Wing sweep angle');
            obj.wingDihedral  = SIM.parameter('Unit','deg','Description','Wing dihedral angle');
            obj.wingIncidence = SIM.parameter('Unit','deg','Description','Wing flow incidence angle');
            obj.wingNACA      = SIM.parameter('Description','Wing NACA airfoil');
            obj.wingClMax     = SIM.parameter('Description','Wing airfoil maximum lift coefficient');
            obj.wingClMin     = SIM.parameter('Description','Wing airfoil minimum lift coefficient');
            % H-stab
            obj.RhsLE_wingLE  = SIM.parameter('Unit','m','Description','Vector going from wing leading edge to H-stab leading edge');
            obj.hsChord     = SIM.parameter('Unit','m','Description','H-stab root chord');
            obj.hsAR        = SIM.parameter('Description','H-stab Aspect ratio');
            obj.hsTR        = SIM.parameter('Description','H-stab Taper ratio');
            obj.hsSweep     = SIM.parameter('Unit','deg','Description','H-stab sweep angle');
            obj.hsDihedral  = SIM.parameter('Unit','deg','Description','H-stab dihedral angle');
            obj.hsIncidence = SIM.parameter('Unit','deg','Description','H-stab flow incidence angle');
            obj.hsNACA      = SIM.parameter('Description','H-stab NACA airfoil');
            obj.hsClMax     = SIM.parameter('Description','H-stab airfoil maximum lift coefficient');
            obj.hsClMin     = SIM.parameter('Description','H-stab airfoil minimum lift coefficient');
            % V-stab
            obj.Rvs_wingLE    = SIM.parameter('Unit','m','Description','Vector going from wing leading edge to V-stab leading edge');
            obj.vsChord     = SIM.parameter('Unit','m','Description','V-stab root chord');
            obj.vsSpan      = SIM.parameter('Unit','m','Description','V-stab span');
            obj.vsTR        = SIM.parameter('Description','V-stab Taper ratio');
            obj.vsSweep     = SIM.parameter('Unit','deg','Description','V-stab sweep angle');
            obj.vsNACA      = SIM.parameter('Description','V-stab NACA airfoil');
            obj.vsClMax     = SIM.parameter('Description','V-stab airfoil maximum lift coefficient');
            obj.vsClMin     = SIM.parameter('Description','V-stab airfoil minimum lift coefficient');
            % initial conditions
            obj.init_inertialCmPos = SIM.parameter('Unit','m','Description','Initial CM position represented in the inertial frame');
            obj.init_inertialCmVel = SIM.parameter('Unit','m/s','Description','Initial CM velocity represented in the inertial frame');
            obj.init_euler         = SIM.parameter('Unit','rad','Description','Initial Euler angles');
            obj.init_angVel        = SIM.parameter('Unit','rad/s','Description','Initial angular velocities');
        end
        
        %% setters
        function setNumTethers(obj,val,units)
            obj.numTethers.setValue(val,units);
        end
        
        function setNumTurbines(obj,val,units)
            obj.numTurbines.setValue(val,units);
        end
        
        function setBuoyFactor(obj,val,units)
            obj.buoyFactor.setValue(val,units);
        end
        
        function setFluidDensity(obj,val,units)
            obj.fluidDensity.setValue(val,units);
        end
        
        function setVolume(obj,val,units)
            obj.volume.setValue(val,units);
        end
        
        function setIxx(obj,val,units)
            obj.Ixx.setValue(val,units);
        end
        
        function setIyy(obj,val,units)
            obj.Iyy.setValue(val,units);
        end
        
        function setIzz(obj,val,units)
            obj.Izz.setValue(val,units);
        end
        
        function setIxy(obj,val,units)
            obj.Ixy.setValue(val,units);
        end
        
        function setIxz(obj,val,units)
            obj.Ixz.setValue(val,units);
        end
        
        function setIyz(obj,val,units)
            obj.Iyz.setValue(val,units);
        end
        
        function setRcb_cm(obj,val,units)
            obj.Rcb_cm.setValue(reshape(val,3,1),units);
        end
        
        function setRbridle_cm(obj,val,units)
            obj.Rbridle_cm.setValue(reshape(val,3,1),units);
        end
        
        function setFluidCoeffsFileName(obj,val,units)
            if ~endsWith(val,'.mat')
                val = [val '.mat'] ;
            end
           obj.fluidCoeffsFileName.setValue(val,units); 
        end
        
        % wing
        function setRwingLE_cm(obj,val,units)
            obj.RwingLE_cm.setValue(reshape(val,3,1),units);
        end
        
        function setWingChord(obj,val,units)
            obj.wingChord.setValue(val,units);
        end
        
        function setWingAR(obj,val,units)
            obj.wingAR.setValue(val,units);
        end
        
        function setWingTR(obj,val,units)
            obj.wingTR.setValue(val,units);
        end
        
        function setWingSweep(obj,val,units)
            obj.wingSweep.setValue(val,units);
        end
        
        function setWingDihedral(obj,val,units)
            obj.wingDihedral.setValue(val,units);
        end
        
        function setWingIncidence(obj,val,units)
            obj.wingIncidence.setValue(val,units);
        end
        
        function setWingNACA(obj,val,units)
            obj.wingNACA.setValue(val,units);
        end
        
        function setWingClMax(obj,val,units)
            obj.wingClMax.setValue(val,units);
        end
        
        function setWingClMin(obj,val,units)
            obj.wingClMin.setValue(val,units);
        end
        % H-stab
        function setRhsLE_wingLE(obj,val,units)
            obj.RhsLE_wingLE.setValue(reshape(val,3,1),units);
        end

        function setHsChord(obj,val,units)
            obj.hsChord.setValue(val,units);
        end
        
        function setHsAR(obj,val,units)
            obj.hsAR.setValue(val,units);
        end
        
        function setHsTR(obj,val,units)
            obj.hsTR.setValue(val,units);
        end
        
        function setHsSweep(obj,val,units)
            obj.hsSweep.setValue(val,units);
        end
        
        function setHsDihedral(obj,val,units)
            obj.hsDihedral.setValue(val,units);
        end
        
        function setHsIncidence(obj,val,units)
            obj.hsIncidence.setValue(val,units);
        end        
        
        function setHsNACA(obj,val,units)
            obj.hsNACA.setValue(val,units);
        end
        
        function setHsClMaxl(obj,val,units)
            obj.hsClMax.setValue(val,units);
        end
        
        function setHsClMin(obj,val,units)
            obj.hsClMin.setValue(val,units);
        end             
        
        % V-stab
        function setRvs_wingLE(obj,val,units)
            obj.Rvs_wingLE.setValue(reshape(val,3,1),units);
        end
        
        function setVsChord(obj,val,units)
            obj.vsChord.setValue(val,units);
        end
        
        function setVsSpan(obj,val,units)
            obj.vsSpan.setValue(val,units);
        end             
        
        function setVsTR(obj,val,units)
            obj.vsTR.setValue(val,units);
        end
        
        function setVsSweep(obj,val,units)
            obj.vsSweep.setValue(val,units);
        end
        
        function setVsNACA(obj,val,units)
            obj.vsNACA.setValue(val,units);
        end             
        
        function setVsClMax(obj,val,units)
            obj.vsClMax.setValue(val,units);
        end
        
        function setVsClMin(obj,val,units)
            obj.vsClMin.setValue(val,units);
        end
        
        % initial conditions
        function setInitialCmPos(obj,val,units)
            obj.init_inertialCmPos.setValue(reshape(val,3,1),units);
        end             

        function setInitialCmVel(obj,val,units)
            obj.init_inertialCmVel.setValue(reshape(val,3,1),units);
        end
        
        function setInitialEuler(obj,val,units)
            obj.init_euler.setValue(reshape(val,3,1),units);
        end
        
        function setInitialAngVel(obj,val,units)
            obj.init_angVel.setValue(reshape(val,3,1),units);
        end
        
        %% getters
        % mass
        function val = get.mass(obj)
            val = SIM.parameter('Value',obj.volume.Value*obj.fluidDensity.Value/...
                obj.buoyFactor.Value,...
                'Unit','kg','Description','Vehicle mass');
        end
        
        % MI
        function val = get.MI(obj)
            val = SIM.parameter('Value',[obj.Ixx.Value -abs(obj.Ixy.Value) -abs(obj.Ixz.Value);...
                -abs(obj.Ixy.Value) obj.Iyy.Value -abs(obj.Iyz.Value);...
                -abs(obj.Ixz.Value) -abs(obj.Iyz.Value) obj.Izz.Value],'Unit','kg*m^2',....
                'Description','Moment of inertia matrix');
        end
        
        % added mass
        function val = get.addedMass(obj)
            % dummy variables
            density = obj.fluidDensity.Value;
            chord = obj.wingChord.Value;
            span = chord*obj.wingAR.Value;
            HS_chord = obj.hsChord.Value;
            HS_span = HS_chord*obj.hsAR.Value;
            VS_chord = obj.vsChord.Value;
            VS_span = obj.vsSpan.Value;
            
            % calculate
            m_added_x = pi*density*(span*(0.15*chord/2)^2 + ...
                HS_span*(0.15*HS_chord/2)^2 + VS_span*(0.15*VS_chord/2)^2);
            m_added_y = pi*density*(1.98*span*(chord/2)^2 + ...
                1.98*HS_span*(HS_chord/2)^2 + VS_span*(VS_chord/2)^2);
            m_added_z = pi*density*(span*(chord/2)^2 + ...
                HS_span*(HS_chord/2)^2 + 1.98*VS_span*(VS_chord/2)^2);
            
            % store
            val = SIM.parameter('Value',[m_added_x 0 0;0 m_added_y 0; 0 0 m_added_z],...
                'Unit','kg','Description','Added mass of the system in the body frame');
            
        end
        
        % surface outlines
        function val = get.surfaceOutlines(obj)
            % dummy variables
            R_wle = obj.RwingLE_cm.Value;
            
            w_cr = obj.wingChord.Value;
            w_s = w_cr*obj.wingAR.Value;
            w_ct = w_cr*obj.wingTR.Value;
            w_sweep = obj.wingSweep.Value;
            w_di = obj.wingDihedral.Value;
            
            R_hsle = obj.RhsLE_wingLE.Value;
            hs_cr = obj.hsChord.Value;
            hs_s = hs_cr*obj.hsAR.Value;
            hs_ct = hs_cr*obj.hsTR.Value;
            hs_sweep = obj.hsSweep.Value;
            hs_di = obj.hsDihedral.Value;
            
            R_vsle = obj.Rvs_wingLE.Value;
            vs_cr = obj.vsChord.Value;
            vs_s = obj.vsSpan.Value;
            vs_ct = vs_cr*obj.vsTR.Value;
            vs_sweep = obj.vsSweep.Value;
            
            port_wing =  repmat(R_wle',5,1) +  [0, 0, 0;...
                w_s*tand(w_sweep)/2, -w_s/2, tand(w_di)*w_s/2;...
                (w_s*tand(w_sweep)/2)+w_ct, -w_s/2, tand(w_di)*w_s/2;...
                w_cr, 0, 0;...
                0, 0, 0];
            
            stbd_wing = port_wing.*[ones(5,1),-1*ones(5,1),ones(5,1)];
            
            port_hs = repmat(R_wle',5,1) + repmat(R_hsle',5,1) + [0, 0, 0;...
                hs_s*tand(hs_sweep)/2, -hs_s/2, tand(hs_di)*hs_s/2;...
                (hs_s*tand(hs_sweep)/2)+hs_ct,   -hs_s/2, 0;...
                hs_cr, 0, 0;...
                0, 0, 0];
            
            stbd_hs = port_hs.*[ones(5,1),-1*ones(5,1),ones(5,1)];
            
            top_vs = repmat(R_wle',5,1) + repmat(R_vsle',5,1) + [0, 0, 0;...
                vs_s*tand(vs_sweep), 0, vs_s;...
                (vs_s*tand(vs_sweep))+vs_ct, 0, vs_s;...
                vs_cr, 0, 0;...
                0, 0, 0];
            
            fuselage = [R_wle';(R_wle+R_vsle)'];
            
            val.port_wing = SIM.parameter('Value',port_wing','Unit','m',...
                'Description','Port wing surface co-ordinates');
            
            val.stbd_wing = SIM.parameter('Value',stbd_wing','Unit','m',...
                'Description','Starboard wing surface co-ordinates');
            
            val.port_hs = SIM.parameter('Value',port_hs','Unit','m',...
                'Description','Port H-stab surface co-ordinates');
            
            val.stbd_hs = SIM.parameter('Value',stbd_hs','Unit','m',...
                'Description','Starboard H-stab surface co-ordinates');
            
            val.top_vs = SIM.parameter('Value',top_vs','Unit','m',...
                'Description','V-stab surface co-ordinates');
            
            val.fuselage = SIM.parameter('Value',fuselage','Unit','m',...
                'Description','Fuselage line co-ordinates');
            
        end
        
        % Tether attachment points
        function val = get.thrAttchPts(obj)
            
           switch obj.numTethers.Value
               case 1
                   val = SIM.parameter('Value',obj.Rbridle_cm.Value,...
                       'Unit','m','Description','Vehicle tether attachment point');
               case 3
                   port_thr = obj.surfaceOutlines.port_wing.Value(:,2);
%                        + [obj.wingChord.Value*obj.wingTR.Value/2;0;0];
                   
                   stbd_thr = port_thr.*[1;-1;1];
                   
                   aft_thr = obj.RwingLE_cm.Value + ...
                       [max(obj.RhsLE_wingLE.Value(1),obj.Rvs_wingLE.Value(1));0;0]...
                       + [max(obj.hsChord.Value,obj.vsChord.Value);0;0];
                   
                   val = SIM.parameter('Value',[port_thr,aft_thr,stbd_thr],...
                       'Unit','m','Description','Vehicle tether attachment point');
           end
           
        end
        
        % turbine attachment points
        function val = get.turbineAttchPts(obj)
            switch obj.numTurbines.Value
               case 2
                   port_wing = obj.surfaceOutlines.port_wing.Value(:,2);
                   stbd_wing = port_wing.*[1;-1;1];
                   
                   val = SIM.parameter('Value',[port_wing,stbd_wing],...
                       'Unit','m','Description','Vehicle turbine attachment point');
                otherwise
                    fprintf('get method not programmed for %d turbines',obj.numTurbines.Value)
                   
           end
            
        end
        
        % aerodynamic forces moment arms
        function val = get.fluidMomentArms(obj)
            portWingArm = obj.surfaceOutlines.port_wing.Value(:,2).*[0;0.5;0.5] +...
                obj.RwingLE_cm.Value + [obj.wingChord.Value*obj.wingAR.Value*tand(obj.wingSweep.Value)/4;0;0] + ...
                [obj.wingChord.Value*(obj.wingTR.Value+1)/8;0;0];
            
            stbdWingArm = portWingArm.*[1;-1;1];
            
            hsArm = obj.surfaceOutlines.port_hs.Value(:,1) + ...
                [obj.hsChord.Value/4;0;0];
            
            vsArm = obj.surfaceOutlines.top_vs.Value(:,2).*[0;0;0.5] + ...
                obj.surfaceOutlines.top_vs.Value(:,1) + [obj.vsSpan.Value*tand(obj.vsSweep.Value)/2;0;0] + ...
                [obj.vsChord.Value*(obj.vsTR.Value+1)/8;0;0];
            
            
            val = SIM.parameter('Value',[portWingArm,stbdWingArm,hsArm,vsArm],'Unit','m',...
                'Description','Fluid dynamic surface moment arms');
            
            
        end
        
        % aerodynamic reference area
        function val = get.fluidRefArea(obj)
            Sref = obj.wingAR.Value*obj.wingChord.Value^2;
            
            val = SIM.parameter('Value',Sref,'Unit','m^2',...
                'Description','Reference area for aerodynamic calculations');
        end
        
        %% other methods
        
        % scale vehcile
        function obj = scale(obj,lengthScaleFactor,densityScaleFactor)
            
            props = findAttrValue(obj,'SetAccess','private');
            for ii = 1:numel(props)
                obj.(props{ii}).scale(lengthScaleFactor,densityScaleFactor);
            end
        end
        
        % fluid dynamic coefficient data
        function calcFluidDynamicCoefffs(obj)
            
            fileLoc = which(obj.fluidCoeffsFileName.Value);
            presFolder = pwd;
            
            
            
            if isfile(fileLoc)
            load(fileLoc,'aeroStruct');
            obj.fluidCoeffData = aeroStruct;
            
            
            else
                fprintf([' The file containing the fluid dynamic coefficient data file does not exist.\n',...
                    ' Would you like to run AVL and create data file ''%s'' ?\n'],obj.fluidCoeffsFileName.Value);
                
                str = input('(Y/N): \n','s');
                
                if isempty(str)
                    str = 'y';
                end
               
                if strcmpi(str,'Y')
                    avlCreateInputFilePart(obj)
                    
                    %% wing
                    alp_max = 55;
                    alp_min = -55;
                    n_steps = 71;
                    % set run cases
                    alphas   = linspace(alp_min,alp_max,n_steps);
                    ailerons = 0;
                    
                    % run AVL for right wing
                    avlProcessPart(obj,'wing',alphas,ailerons,'Parallel',true);
                    load('resultFile','results');
                    
                    [CLWingTab,CDWingTab] = avlPartitionedLookupTable(results);
                    
                    % get wing aileron gains
                    n_case = 10;
                    alphas = 0;
                    ailerons = linspace(-5,5,n_case);
                    
                    avlProcessPart(obj,'wing',alphas,ailerons,'Parallel',true);
                    load('resultFile','results');
                    
                    CL_w = NaN(1,n_case);
                    CD_w = NaN(1,n_case);
                    for ii = 1:n_case
                        CL_w(ii) = results{1}(ii).FT.CLtot;
                        CD_w(ii) = results{1}(ii).FT.CDtot;
                    end
                    
                    CL_kWing = polyfit(ailerons,CL_w,2);
                    CL_kWing(end) = 0;
                    CD_kWing = polyfit(ailerons,CD_w,2);
                    CD_kWing(end) = 0;
                    
                    % port wing data
                    aeroStruct(1).CL = reshape(CLWingTab.Table.Value,[],1);
                    aeroStruct(1).CD = reshape(CDWingTab.Table.Value,[],1);
                    aeroStruct(1).alpha = reshape(CDWingTab.Breakpoints.Value,[],1);
                    aeroStruct(1).GainCL = reshape(CL_kWing,1,[]);
                    aeroStruct(1).GainCD =  reshape(CD_kWing,1,[]);
                    aeroStruct(1).sCb = eye(3);
                    
                    % stbd wing data
                    aeroStruct(2).CL = aeroStruct(1).CL;
                    aeroStruct(2).CD = aeroStruct(1).CD;
                    aeroStruct(2).alpha = aeroStruct(1).alpha;
                    aeroStruct(2).GainCL = aeroStruct(1).GainCL;
                    aeroStruct(2).GainCD =  aeroStruct(1).GainCD;
                    aeroStruct(2).sCb = aeroStruct(1).sCb;
                    
                    %% horizontal stabilizers
                    % set run cases
                    alphas   = linspace(alp_min,alp_max,n_steps);
                    ailerons = 0;
                    
                    % run AVL for HS
                    avlProcessPart(obj,'H_stab',alphas,ailerons,'Parallel',true);
                    load('resultFile','results');
                    
                    [CLHSTab,CDHSTab] = avlPartitionedLookupTable(results);
                    
                    % get HS aileron gains
                    alphas = 0;
                    ailerons = linspace(-5,5,n_case);
                    
                    avlProcessPart(obj,'H_stab',alphas,ailerons,'Parallel',true);
                    load('resultFile','results');
                    
                    CL_hs = NaN(1,n_case);
                    CD_hs = NaN(1,n_case);
                    
                    for ii = 1:n_case
                        CL_hs(ii) = results{1}(ii).FT.CLtot;
                        CD_hs(ii) = results{1}(ii).FT.CDtot;
                    end
                    
                    CL_kHS = polyfit(ailerons,CL_hs,2);
                    CL_kHS(end) = 0;
                    CD_kHS = polyfit(ailerons,CD_hs,2);
                    CD_kHS(end) = 0;
                    
                    % HS data
                    aeroStruct(3).CL = reshape(CLHSTab.Table.Value,[],1);
                    aeroStruct(3).CD = reshape(CDHSTab.Table.Value,[],1);
                    aeroStruct(3).alpha = reshape(CDHSTab.Breakpoints.Value,[],1);
                    aeroStruct(3).GainCL = reshape(CL_kHS,1,[]);
                    aeroStruct(3).GainCD =  reshape(CD_kHS,1,[]);
                    aeroStruct(3).sCb = eye(3);
                    
                    %% vertical stabilizer
                    % set run cases
                    alphas   = linspace(alp_min,alp_max,n_steps);
                    ailerons = 0;
                    
                    % run AVL for VS
                    avlProcessPart(obj,'V_stab',alphas,ailerons,'Parallel',true);
                    load('resultFile','results');
                    
                    [CLVSTab,CDVSTab] = avlPartitionedLookupTable(results);
                    
                    % get VS aileron gains
                    alphas = 0;
                    ailerons = linspace(-5,5,n_case);
                    
                    avlProcessPart(obj,'V_stab',alphas,ailerons,'Parallel',true);
                    load('resultFile','results');
                    
                    CL_vs = NaN(1,n_case);
                    CD_vs = NaN(1,n_case);
                    for ii = 1:n_case
                        CL_vs(ii) = results{1}(ii).FT.CLtot;
                        CD_vs(ii) = results{1}(ii).FT.CDtot;
                    end
                    
                    CL_kVS = polyfit(ailerons,CL_vs,2);
                    CL_kVS(end) = 0;
                    CD_kVS = polyfit(ailerons,CD_vs,2);
                    CD_kVS(end) = 0;
                    
                    aeroStruct(4).CL = reshape(CLVSTab.Table.Value,[],1);
                    aeroStruct(4).CD = reshape(CDVSTab.Table.Value,[],1);
                    aeroStruct(4).alpha = reshape(CDVSTab.Breakpoints.Value,[],1);
                    aeroStruct(4).GainCL = reshape(CL_kVS,1,[]);
                    aeroStruct(4).GainCD =  reshape(CD_kVS,1,[]);
                    aeroStruct(4).sCb = [1 0 0;0 0 -1;0 1 0];
                    
                    aeroStruct = reshape(aeroStruct,1,[]);
                    
                    obj.fluidCoeffData = aeroStruct;
                    
                    save(obj.fluidCoeffsFileName.Value,'aeroStruct');

                    delete('wing');
                    delete('H_stab');
                    delete('V_stab');
                    
                    %
                    filepath = fileparts(which('avl.exe'));
                    
                    delete(fullfile(filepath,strcat('resultFile','.mat')));
                    
                    fprintf('''%s'' created in:\n %s\n',...
                        obj.fluidCoeffsFileName.Value,fileparts(which(obj.fluidCoeffsFileName.Value)));
                    
                    cd(presFolder);
                    
                elseif strcmpi(str,'N')
                    warning('Simulation won''t run without valid aero coefficient values')
                    
                else
                    error('Invalid input')
                end

            end
  
        end
        
        % determine aerodynamic center of the design
        function calcDesginFluidDynamicCenter(obj)
            [~,ii_max] = max(obj.fluidCoeffData(1).CL);
            [~,ii_min] = min(abs(obj.fluidCoeffData(1).CL));
            alphas = obj.fluidCoeffData(1).alpha(ii_min:ii_max);
            
            for ii = 1:3
                CLs(:,ii) = obj.fluidCoeffData(ii).CL(ii_min:ii_max,1);
                CDs(:,ii) = obj.fluidCoeffData(ii).CD(ii_min:ii_max,1);
                
            end
            
            init_RwingLE_cm = obj.RwingLE_cm.Value;
            
            base_portWingArm = obj.fluidMomentArms.Value(:,1)-init_RwingLE_cm;
            base_stbdWingArm = obj.fluidMomentArms.Value(:,2)-init_RwingLE_cm;
            base_HSArm = obj.fluidMomentArms.Value(:,3)-init_RwingLE_cm;
            
            baseArms = [base_portWingArm,base_stbdWingArm,base_HSArm];
            
            w_span = obj.wingAR.Value*obj.wingChord.Value;
            
                        
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            lb = 0;
            ub = max( max(obj.surfaceOutlines.port_hs.Value(1,:)),...
                max(obj.surfaceOutlines.top_vs.Value(1,:)) );
            
            
            [x_fin,fval] = fmincon(...
                @(x) calcCM(x,baseArms,w_span,CLs,CDs,alphas,1),init_RwingLE_cm(1),...
                A,b,Aeq,beq,lb,ub);
            
            
            CM_all = calcCM(x_fin,baseArms,w_span,CLs,CDs,alphas,0);
            
        end
        
        % plotting functions
        function plot(obj)
            
            fh = findobj( 'Type', 'Figure', 'Name', 'Design');
            
            if isempty(fh)
                fh = figure;
                fh.Name ='Design';
            else
                figure(fh);
            end
            
            fs = fieldnames(obj.surfaceOutlines);
            
            for ii = 1:6
                plot3(obj.surfaceOutlines.(fs{ii}).Value(1,:),...
                    obj.surfaceOutlines.(fs{ii}).Value(2,:),...
                    obj.surfaceOutlines.(fs{ii}).Value(3,:),...
                    'LineWidth',1.2,'Color','k','LineStyle','-');
                hold on
            end
            
            for ii = 1:obj.numTethers.Value
                pTet = plot3(obj.thrAttchPts.Value(1,ii),...
                    obj.thrAttchPts.Value(2,ii),...
                    obj.thrAttchPts.Value(3,ii),...
                    'r+');
                
            end
            
            for ii = 1:obj.numTurbines.Value
                pTurb = plot3(obj.turbineAttchPts.Value(1,ii),...
                    obj.turbineAttchPts.Value(2,ii),...
                    obj.turbineAttchPts.Value(3,ii),...
                    'm+');
            end
            
            for ii = 1:4
                pMom = plot3(obj.fluidMomentArms.Value(1,ii),...
                    obj.fluidMomentArms.Value(2,ii),...
                    obj.fluidMomentArms.Value(3,ii),...
                    'b+');
                
            end
            
            pCM = plot3(0,0,0,'r*');
            grid on
            axis equal
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            view(-45,30)
            legend([pCM,pTet,pTurb,pMom],{'CM','Tethered pts.','Turbine pts.','Aero force pts.'},...
                'Location','northeast')
            
        end
        
        function plotCoeffPolars(obj)
            fh = findobj( 'Type', 'Figure', 'Name', 'Partitioned Aero Coeffs');
            
            if isempty(fh)
                fh = figure;
                fh.Position =[102 92 3*560 2*420];
                fh.Name ='Partitioned Aero Coeffs';
            else
                figure(fh);
            end
            
            % left wing
            ax1 = subplot(2,4,1);
            plot(obj.fluidCoeffData(1).alpha,obj.fluidCoeffData(1).CL);
            hCL_ax = gca;
            
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('Port Wing')
            grid on
            hold on
            
            ax5 = subplot(2,4,5);
            plot(obj.fluidCoeffData(1).alpha,obj.fluidCoeffData(1).CD);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            hCD_ax = gca;
            
            linkaxes([ax1,ax5],'x');
            
            % right wing
            ax2 = subplot(2,4,2);
            plot(obj.fluidCoeffData(2).alpha,obj.fluidCoeffData(2).CL);
            
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('Stbd Wing')
            grid on
            hold on
            
            ax6 = subplot(2,4,6);
            plot(obj.fluidCoeffData(2).alpha,obj.fluidCoeffData(2).CD);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            
            linkaxes([ax2,ax6],'x');
            
            % HS
            ax3 = subplot(2,4,3);
            plot(obj.fluidCoeffData(3).alpha,obj.fluidCoeffData(3).CL);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('H-stab')
            grid on
            hold on
            
            ax7 = subplot(2,4,7);
            plot(obj.fluidCoeffData(3).alpha,obj.fluidCoeffData(3).CD);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            
            linkaxes([ax3,ax7],'x');
            
            % VS
            ax4 = subplot(2,4,4);
            plot(obj.fluidCoeffData(4).alpha,obj.fluidCoeffData(4).CL);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{L}$')
            title('V-stab')
            grid on
            hold on
            
            ax8 = subplot(2,4,8);
            plot(obj.fluidCoeffData(4).alpha,obj.fluidCoeffData(4).CD);
            xlabel('$\alpha$ [deg]')
            ylabel('$C_{D}$')
            grid on
            hold on
            
            linkaxes([ax4,ax8],'x');
            
            axis([ax1 ax2 ax3 ax4],[-inf inf hCL_ax.YLim(1) hCL_ax.YLim(2)]);
            axis([ax5 ax6 ax7 ax8],[-inf inf hCD_ax.YLim(1) hCD_ax.YLim(2)]);
            
        end
        

        
    end
    
    
    
    
end
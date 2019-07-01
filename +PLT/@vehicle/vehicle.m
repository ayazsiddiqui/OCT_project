classdef vehicle < dynamicprops
    
    %VEHICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        lengthScale
        densityScale
        numTethers
        numTurbines
        buoyFactor
        volume
        Ixx
        Iyy
        Izz
        Ixy
        Ixz
        Iyz
        Rcb_cm
        % aero properties
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
    
    properties (Dependent)
        mass
        
    end
    
    methods
        %% Constructor
        function obj = vehicle
            %VEHICLE Construct an instance of this class
            obj.lengthScale  = SIM.parameter('Description','Length scale factor');
            obj.densityScale = SIM.parameter('Description','Length scale factor');
            obj.numTethers  = SIM.parameter('Description','Number of tethers');
            obj.numTurbines = SIM.parameter('Description','Number of turbines');
            obj.buoyFactor = SIM.parameter('Description','Buoyancy Factor');
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
        function setLengthScale(obj,val,units)
            obj.lengthScale.setValue(val,units);
        end
        
        function setDensityScale(obj,val,units)
            obj.densityScale.setValue(val,units);
        end
        
        function setNumTethers(obj,val,units)
            obj.numTethers.setValue(val,units);
        end
        
        function setNumTurbines(obj,val,units)
            obj.numTurbines.setValue(val,units);
        end
        
        function setBuoyFactor(obj,val,units)
            obj.buoyFactor.setValue(val,units);
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
            obj.Rcb_cm.setValue(val,units);
        end
        
        % wing
        function setRwingLE_cm(obj,val,units)
            obj.RwingLE_cm.setValue(val,units);
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
            obj.RhsLE_wingLE.setValue(val,units);
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
            obj.Rvs_wingLE.setValue(val,units);
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
            obj.init_inertialCmPos.setValue(val,units);
        end             

        function setInitialCmVel(obj,val,units)
            obj.init_inertialCmVel.setValue(val,units);
        end
        
        function setInitialEuler(obj,val,units)
            obj.init_euler.setValue(val,units);
        end
        
        function setInitialAngVel(obj,val,units)
            obj.init_angVel.setValue(val,units);
        end
        
        %% getters
        function val = get.mass(obj)
            val = SIM.parameter('Value',1e3*obj.volume.Value*obj.densityScale.Value/...
                obj.buoyFactor.Value,...
                'Unit','kg','Description','Vehicle mass');
        end
        
        
        
        %% other methods
        function plot(obj)
            
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
            
            all_coords = {port_wing;stbd_wing;port_hs;stbd_hs;top_vs;fuselage};
            
            for ii = 1:6
                plot3(all_coords{ii}(:,1),...
                    all_coords{ii}(:,2),...
                    all_coords{ii}(:,3),...
                    'LineWidth',1.2,'Color','k','LineStyle','-')
                hold on
            end
            plot3(0,0,0,'r*')
            grid on
            axis equal
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
          
            
            
        end
        
        function scaleVehicle(obj)
            LS = obj.lengthScale.Value;
            
            % scale volume and inetias
            obj.setVolume(obj.volume.Value*LS^3,'m^3');
            obj.setIxx(obj.Ixx.Value*LS^5,'kg*m^2');
            obj.setIyy(obj.Iyy.Value*LS^5,'kg*m^2');
            obj.setIzz(obj.Izz.Value*LS^5,'kg*m^2');
            obj.setIxy(obj.Ixy.Value*LS^5,'kg*m^2');
            obj.setIxz(obj.Ixz.Value*LS^5,'kg*m^2');
            obj.setIyz(obj.Iyz.Value*LS^5,'kg*m^2');
            obj.setRcb_cm(obj.Rcb_cm.Value*LS,'m');
            
            % scale wing
            obj.setRwingLE_cm(obj.RwingLE_cm.Value*LS,'m');
            obj.setWingChord(obj.wingChord.Value*LS,'m');
            
            % scale H-stab
            obj.setRhsLE_wingLE(obj.RhsLE_wingLE.Value*LS,'m');
            obj.setHsChord(obj.hsChord.Value*LS,'m');
            
            % sacle V-stab
            obj.setRvs_wingLE(obj.Rvs_wingLE.Value*LS,'m');
            obj.setVsChord(obj.vsChord.Value*LS,'m');
            obj.setVsSpan(obj.vsSpan.Value*LS,'m');
            obj.setVsTR(0.8,'');
            
            % initial conditions
            obj.setInitialCmPos(obj.init_inertialCmPos.Value.*LS,'m');
            obj.setInitialCmVel(obj.init_inertialCmVel.Value.*LS^0.5,'m/s');
            obj.setInitialAngVel(obj.init_angVel.Value.*(1/LS^0.5),'rad/s');
            
        end
            
            
            
            
            
            
            
        
        
        

        
    end
    
    
    
    
end
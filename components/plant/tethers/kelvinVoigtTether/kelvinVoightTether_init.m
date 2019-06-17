% Kelvin Voigt tether initialization scrit
numTethers = length(thr_struct);

% initialize variables
numNodes = zeros(1,numTethers);
diameters = zeros(1,numTethers);
youngsModuli = zeros(1,numTethers);
dampingRatios = zeros(1,numTethers);
dragCoeffs = zeros(1,numTethers);
densities = zeros(1,numTethers);
vehicleMasses = zeros(1,numTethers);

for ii = 1:numTethers
    
    numNodes(ii) = thr_struct(ii).numNodes;
    diameters(ii) = thr_struct(ii).diameter;
    youngsModuli(ii) = thr_struct(ii).youngsModulus;
    dampingRatios(ii) = thr_struct(ii).dampingRatio;
    dragCoeffs(ii) = thr_struct(ii).dragCoeff;
    densities(ii) = thr_struct(ii).density;
    vehicleMasses(ii) = thr_struct(ii).vehicleMass;
    
end

numNodes = min(numNodes);
youngsModuli = min(youngsModuli);
dampingRatios = min(dampingRatios);
dragCoeffs = min(dragCoeffs);
densities = min(densities);
vehicleMasses = min(vehicleMasses);





% The active variant can only be modified at edit-time or very early during
% the simulink compile stage.
% simulationStatus =  get_param(bdroot,'SimulationStatus');
% 
% vssBlk = gcb;
% oldVariant = get_param(gcb,'OverrideUsingVariant');
% 
% if any(numNodes==2)
%     set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNode2_variant');
% else
%     set_param(vssBlk, 'OverrideUsingVariant', 'tetherNumNodeN_variant');
% end

    
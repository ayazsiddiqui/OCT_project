function bayesianAscentSfunc(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 3;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties

% design variables
block.InputPort(1).Dimensions        = [size(block.DialogPrm(1).Data,1),1];
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;

% design fval
block.InputPort(2).Dimensions        = [1,1];
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;

% iteration number
block.InputPort(3).Dimensions        = [1,1];
block.InputPort(3).DatatypeID  = 0;  % double
block.InputPort(3).Complexity  = 'Real';
block.InputPort(3).DirectFeedthrough = true;

% Override output port properties
block.OutputPort(1).Dimensions       = block.InputPort(1).Dimensions;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';


% Register parameters
block.NumDialogPrms     = 6;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)

trainDsgn =  block.DialogPrm(1).Data;
noInputs = size(trainDsgn,1);
trainFval =  block.DialogPrm(2).Data;
trainOpHyp =  block.DialogPrm(3).Data;
noiseVar = block.DialogPrm(4).Data;
designLimits = block.DialogPrm(5).Data;
maxIter = block.DialogPrm(6).Data;

iniTau = 1.1;
gamma = 0.01;
beta = 1.1;

x0 = block.InputPort(1).Data;
xFval = block.InputPort(2).Data;
noIter = block.InputPort(3).Data;

% preallocate matrices
noTrainDsgn = numel(trainFval);
nt = noTrainDsgn;
testDsgns = NaN(noInputs,nt+maxIter);
testFval = NaN(nt + maxIter,1);
testOpHyp = NaN(size(trainOpHyp,1),maxIter);
finPts = NaN(noInputs,maxIter);
finFval = NaN(maxIter,1);
tau = NaN(size(iniTau,1),maxIter);

% intitialize
testDsgns(:,1:nt) = trainDsgn;
testFval(1:nt,1) = trainFval;
testOpHyp(:,1) = trainOpHyp;

while noIter <= maxIter
    if noIter == 1
        testDsgns(:,nt+noIter) = x0;
        testFval(nt+noIter,1) = xFval;
        testOpHyp(:,noIter) = trainOpHyp;
        finPts(:,noIter) = x0;
        finFval(noIter,1) = testFval(nt+noIter,1);
        tau(:,noIter) = iniTau;
        
    else
        testDsgns(:,nt+noIter) = optPt;
        testFval(nt+noIter,1) = optFval;
        
        % optimize hyper parameters
        A = []; b = [];
        Aeq = []; beq = [];
        
        % bounds
        lb = [eps,1e-2*ones(block.InputPort(1).Dimensions)'];
        ub = [10,10*ones(block.InputPort(1).Dimensions)'];
        nonlcon = [];
        options  = optimoptions('fmincon','Display','off');
        
        optHyp = fmincon(@(hyper) ...
            -calcLogLikelihood(trainDsgn,trainFval,hyper),...
            trainHyp,A,b,Aeq,beq,lb,ub,nonlcon,options);
        
        testOpHyp(:,noIter) = OpHyp;
        finPts(:,noIter) = optPt;
        finFval(noIter,1) = optFval;
        
        if finFval(noIter,1)-finFval(noIter-1,1) >= gamma*(1/noIter)*(max(testFval(1:noIter-1,1))-finFval(1))
            tau(:,noIter) = beta*tau(:,noIter-1);
            
        else
            tau(:,noIter) = iniTau;
        end
    end
    
    % step 2: construct GP model
    testCovMat = buildCovarianceMatrix(testDsgns(:,1:nt+noIter),testDsgns(:,1:nt+noIter),testOpHyp(1,noIter),noiseVar,testOpHyp(2:end,noIter));
    
    % select next input
    xLims = calDesignBounds(finPts(:,noIter),tau(:,noIter),designLimits);
    
    
    % maximize acquisition function
    [optX,~] = maximizeAcquisitionFunction(max(finFval),testDsgns(:,1:nt+noIter),testCovMat,...
    testFval(1:nt+noIter,1),finPts(:,noIter),xLims,noiseVar,testOpHyp(:,noIter));

end


% opt
block.OutputPort(1).Data = optX;


%end Outputs

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate



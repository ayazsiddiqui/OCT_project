function bayesianAscentSfunc2(block)
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
block.NumInputPorts  = 2;
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

% Override output port properties
block.OutputPort(1).Dimensions       = block.InputPort(1).Dimensions;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';

% Register parameters
block.NumDialogPrms     = 10;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [-1, 0];

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
block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)
block.NumDworks = 7;

trainDsgn =  block.DialogPrm(1).Data;
noInputs = size(trainDsgn,1);
trainFval =  block.DialogPrm(2).Data;
maxIter = ceil(block.DialogPrm(7).Data) + 1;

block.Dwork(1).Name            = 'testDsgns';
block.Dwork(1).Dimensions      = noInputs*(length(trainFval) + maxIter);
block.Dwork(1).DatatypeID      = 0;      % double
block.Dwork(1).Complexity      = 'Real'; % real
block.Dwork(1).UsedAsDiscState = true;

block.Dwork(2).Name            = 'testFval';
block.Dwork(2).Dimensions      = 1*(length(trainFval) + maxIter);
block.Dwork(2).DatatypeID      = 0;      % double
block.Dwork(2).Complexity      = 'Real'; % real
block.Dwork(2).UsedAsDiscState = true;

block.Dwork(3).Name            = 'testOpHyp';
block.Dwork(3).Dimensions      = (noInputs+1)*(maxIter);
block.Dwork(3).DatatypeID      = 0;      % double
block.Dwork(3).Complexity      = 'Real'; % real
block.Dwork(3).UsedAsDiscState = true;

block.Dwork(4).Name            = 'finPts';
block.Dwork(4).Dimensions      = noInputs*(maxIter);
block.Dwork(4).DatatypeID      = 0;      % double
block.Dwork(4).Complexity      = 'Real'; % real
block.Dwork(4).UsedAsDiscState = true;

block.Dwork(5).Name            = 'finFval';
block.Dwork(5).Dimensions      = 1*(maxIter);
block.Dwork(5).DatatypeID      = 0;      % double
block.Dwork(5).Complexity      = 'Real'; % real
block.Dwork(5).UsedAsDiscState = true;

block.Dwork(6).Name            = 'tau';
block.Dwork(6).Dimensions      = noInputs*(maxIter);
block.Dwork(6).DatatypeID      = 0;      % double
block.Dwork(6).Complexity      = 'Real'; % real
block.Dwork(6).UsedAsDiscState = true;

block.Dwork(7).Name            = 'noIter';
block.Dwork(7).Dimensions      = 1;
block.Dwork(7).DatatypeID      = 0;      % double
block.Dwork(7).Complexity      = 'Real'; % real
block.Dwork(7).UsedAsDiscState = true;


%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)
trainDsgn =  block.DialogPrm(1).Data;
noInputs = size(trainDsgn,1);
trainFval =  block.DialogPrm(2).Data;

block.Dwork(1).Data(1:noInputs*length(trainFval)) = trainDsgn(:);
block.Dwork(2).Data(1:1*length(trainFval)) = trainFval(:);
block.Dwork(7).Data = 1;


%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)

noInputs = size(block.DialogPrm(1).Data,1);
nTrain =  size(block.DialogPrm(1).Data,2);
designLimits = block.DialogPrm(3).Data;

gamma = block.DialogPrm(4).Data;
beta = block.DialogPrm(5).Data;
iniTau = block.DialogPrm(6).Data*(max(designLimits,[],2)-min(designLimits,[],2));

x0 = block.InputPort(1).Data;
xFval = block.InputPort(2).Data;
noIter = block.Dwork(7).Data;
gp = block.DialogPrm(8).Data;
predHorizon = block.DialogPrm(9).Data;
ctrlHorizon = block.DialogPrm(10).Data;


block.Dwork(1).Data(noInputs*(nTrain + noIter -1)+1:noInputs*(nTrain + noIter)) = x0;
block.Dwork(2).Data(nTrain + noIter) = xFval;
block.Dwork(4).Data(noInputs*(noIter -1)+1:noInputs*(noIter)) = x0;
block.Dwork(5).Data(noIter) = xFval;

% allocatung dummy variables
trainDsgns = reshape(block.Dwork(1).Data(1:noInputs*nTrain),...
    noInputs,[]);
trainFval = reshape(block.Dwork(2).Data(1:nTrain),[],1);

finPts = reshape(block.Dwork(4).Data(1:noInputs*noIter),...
    noInputs,[]);
finFval = reshape(block.Dwork(5).Data(1:noIter),[],1);

if noIter == 1
    OpHyp = rand(1+noInputs,1);
    tau = iniTau;
else
    OpHyp = reshape(block.Dwork(3).Data(1:(noInputs+1)*noIter),noInputs+1,[]);
    tau = reshape(block.Dwork(6).Data(1:(noInputs*noIter)),noInputs,[]);
end

[sol] = gp.mpcBayesianAscent(trainDsgns,trainFval,finPts,finFval,...
    OpHyp,tau,designLimits,iniTau,gamma,beta,noIter,predHorizon,ctrlHorizon);

block.Dwork(3).Data((noInputs+1)*(noIter-1)+1:(noInputs+1)*(noIter)) = sol.testOpHyp;
block.Dwork(6).Data((noInputs)*(noIter-1)+1:(noInputs*noIter)) = sol.tau(:,noIter);


% opt
block.OutputPort(1).Data = sol.optPt;
block.Dwork(7).Data = block.Dwork(7).Data + 1;


%end Outputs

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

% Terminate(block)

%end Terminate
% gpTrain,gpTrainFval,noiseVars,desLims,gamma,beta,iniTauPerc,maxIter


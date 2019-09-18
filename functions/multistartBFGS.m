function [optDsgn,maxF] = multistartBFGS(objF,X,lb,ub,varargin)
%PARTICLESWARMOPT
% Function to find the global maximum of an objective function using
% Particle Swarm Optimization.
%
% Required inputs:
% objF = objective function which should be a function of X
% X =  dummy design matrix, rows represent inputs and columns represent different
% designs
% lb, ub =  lower and upper bounds on the design, should be the same size
% as X
%
% varargin parameters:
% 'swarmSize' = number of particles exploring the design space
% 'cognitiveLR' = cognitive (individul) learning rate. inidicates the
% partices tendency to gravitate towards the best value it found
% 'socialLR' = social (group) learning rate. indicates the particle's
% tendency to gravitate towards the best value the swarm found
% 'maxIter' = maximum number of iterations
%
% Example use of code
% [maxF,optDsgn] = particleSwarmOpt(@(x)objF(x),X,lb,ub,...
%     'swarmSize',25,'cognitiveLR',0.4,'socialLR',0.2,'maxIter',20);

p = inputParser;
addRequired(p,'objF');
addRequired(p,'X',@isnumeric);
addRequired(p,'lb',@isnumeric);
addRequired(p,'ub',@isnumeric);
addParameter(p,'nStarts',25,@isnumeric);

parse(p,objF,X,lb,ub,varargin{:});

ss = p.Results.swarmSize;

% design space size
dsgnSize = size(X);
% initial swarm size
iniSwarm = lb + (ub-lb).*rand([dsgnSize,ss]);
swarm = NaN(size(iniSwarm));

% initial design space
fVal = NaN(ss,1);

for ii = 1:ss
    [mpcOptPts,mpcOptFval] = BFGS(@(X) -p.Results.objF(X),iniSwarm(:,:,ii),...
        'lb',lb,'ub',ub,'maxIter',20,'bfgsConvergeTol',1e-2,'bpStep',0.25,...
        'bpMaxIter',100,'gradStep',0.02,'GsConvergeTol',5e-3);
%     
%     options  = optimoptions('fmincon','Display','off');
% 
%     [mpcOptPts,mpcOptFval] = fmincon(@(X) -p.Results.objF(X),iniSwarm(:,:,ii),...
%         [],[],[],[],lb,ub,[],options);
    
    
    fVal(ii) = mpcOptFval;
    swarm(:,:,ii) = mpcOptPts;
    
end

[minima,idx] = min(fVal);
GbesLoc = swarm(:,:,idx);

maxF = -minima;
optDsgn = GbesLoc;

end


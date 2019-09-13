function [optDsgn,maxF] = sequentialParticleSwarmOpt(objF,iniPt,nSteps,lb,ub,varargin)
%PARTICLESWARMOPT
% Function to find the global maximum of an objective function using
% Particle Swarm Optimization.
% 
% Required inputs:
% objF = objective function which should be a function of X
% iniPt = starting point
% nSteps = number of steps to take in sequence
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
addRequired(p,'iniPt',@isnumeric);
addRequired(p,'nSteps',@isnumeric);
addRequired(p,'lb',@isnumeric);
addRequired(p,'ub',@isnumeric);
addParameter(p,'swarmSize',25,@isnumeric);
addParameter(p,'cognitiveLR',0.4,@isnumeric);
addParameter(p,'socialLR',0.2,@isnumeric);
addParameter(p,'maxIter',20,@isnumeric);
addParameter(p,'stepPerc',0.2,@isnumeric);

parse(p,objF,iniPt,nSteps,lb,ub,varargin{:});

ss = p.Results.swarmSize;
maxIter = p.Results.maxIter;

dsgnRange = p.Results.ub - p.Results.lb;
tau = p.Results.stepPerc*dsgnRange;
% design space size
dsgnSize = size(iniPt);
% initial swarm size
iniSwarm = NaN(dsgnSize(1),p.Results.nSteps,ss);
for jj = 1:p.Results.nSteps
    if jj == 1
        iniSwarm(:,jj,:) = (p.Results.iniPt-tau) + (2*tau).*rand([dsgnSize,ss]);
    else
        iniSwarm(:,jj,:) = (iniSwarm(:,jj-1,:)-tau) + (2*tau).*rand([dsgnSize,ss]);
    end
end
% initial design space
fVal = NaN(ss,maxIter);

% swarm
swarm = NaN([size(iniSwarm),maxIter]);
V = NaN(size(swarm));
PbestLoc = NaN(size(iniSwarm));

if size(lb,2) ~= size(iniSwarm,2)
    lb = repmat(lb,1,size(iniSwarm,2));
end
if size(ub,2) ~= size(iniSwarm,2)
    ub = repmat(ub,1,size(iniSwarm,2));
end    
    
for jj = 1:maxIter
    if jj == 1
        V(:,:,:,jj) = zeros(size(iniSwarm));
        swarm(:,:,:,jj) = iniSwarm;
    else
        V(:,:,:,jj) = V(:,:,:,jj-1) + p.Results.cognitiveLR*rand*(PbestLoc-swarm(:,:,:,jj-1))...
            + p.Results.socialLR*rand*(GbesLoc-swarm(:,:,:,jj-1));
        swarm(:,:,:,jj) = V(:,:,:,jj) + swarm(:,:,:,jj-1);
        
        for ii = 1:ss
            [row,col] = find(swarm(:,:,ii,jj)<lb);
            swarm(row,col,ii,jj) = lb(row,col);
            [row,col] = find(swarm(:,:,ii,jj)>ub);
            swarm(row,col,ii,jj) = ub(row,col);
        end
    end
    
    for ii = 1:ss
        fVal(ii,jj) = p.Results.objF(swarm(:,:,ii,jj));
        
        [~,PbestIdx] = max(fVal(ii,:));
        PbestLoc(:,:,ii) = swarm(:,:,ii,PbestIdx(1));
        
    end
    maximum = max(max(fVal));
    [x,y]=find(fVal==maximum);
    GbesLoc = swarm(:,:,x(1),y(1));
    
end

maxF = maximum;
optDsgn = GbesLoc;

end


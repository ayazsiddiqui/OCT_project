function [maxF,optDsgn] = particleSwarmOpt(objF,X,lb,ub,varargin)
%PARTICLESWARMOPT Summary of this function goes here
% function to find the global maximum of an objective function using
% particle swarm optimization

%   Detailed explanation goes here
% objF = objective function which should be a function of X
% X =  dummy design matrix, rows represent inputs and columns represent different
% designs
% lb, ub =  lower and upper bounds on the design, should be the same size
% as X
% varargin
% 'swarmSize' = number of particles exploring the design space
% 'cognitiveLR' = cognitive (individul) learning rate
% 'socialLR' = social (group) learning rate
% 'maxIter' = maximum number of iterations

p = inputParser;
addRequired(p,'objF');
addRequired(p,'X',@isnumeric);
addRequired(p,'lb',@isnumeric);
addRequired(p,'ub',@isnumeric);
addParameter(p,'swarmSize',20,@isnumeric);
addParameter(p,'cognitiveLR',1,@isnumeric);
addParameter(p,'socialLR',1,@isnumeric);
addParameter(p,'maxIter',20,@isnumeric);

parse(p,objF,X,lb,ub,varargin{:});

ss = p.Results.swarmSize;
maxIter = p.Results.maxIter;

% design space size
dsgnSize = size(X);
% initial swarm size
iniSwarm = lb + (ub-lb).*rand([dsgnSize,ss]);
iniSwarm = cat(3,X,iniSwarm);

% initial design space
fVal = NaN(ss+1,maxIter);

% swarm
swarm = NaN([size(iniSwarm),maxIter]);
V = NaN(size(swarm));
PbestLoc = NaN(size(iniSwarm));

for jj = 1:maxIter
    if jj == 1
        V(:,:,:,jj) = zeros(size(iniSwarm));
        swarm(:,:,:,jj) = iniSwarm;
    else
        V(:,:,:,jj) = V(:,:,:,jj-1) + p.Results.cognitiveLR*rand*(PbestLoc-swarm(:,:,:,jj-1))...
            + p.Results.socialLR*rand*(GbesLoc-swarm(:,:,:,jj-1));
        swarm(:,:,:,jj) = V(:,:,:,jj) + swarm(:,:,:,jj-1);
        
        [row,col,v] = find(swarm(:,:,:,jj)<lb);
        swarm(row,col,v,jj) = lb(row,col,v);

        [row,col,v] = find(swarm(:,:,:,jj)>ub);
        swarm(row,col,v,jj) = ub(row,col,v);
    end
    
    for ii = 1:ss+1
        fVal(ii,jj) = p.Results.objF(swarm(:,:,ii,jj));
        
        [~,PbestIdx] = max(fVal(ii,:));
        PbestLoc(:,:,ii) = swarm(:,:,ii,PbestIdx);
        
    end
    maximum = max(max(fVal));
    [x,y]=find(fVal==maximum);
    GbesLoc = swarm(:,:,x,y);
    
end



outputArg1 = inputArg1;
outputArg2 = inputArg2;
end


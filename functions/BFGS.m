function [optPts,fMin] = BFGS(objF,iniPt,varargin)
%BFGS Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addRequired(p,'objF');
addRequired(p,'iniPt',@isnumeric);
addParameter(p,'lb',-inf*iniPt,@isnumeric);
addParameter(p,'ub',inf*iniPt,@isnumeric);
addParameter(p,'maxIter',100,@isnumeric);
addParameter(p,'bfgsConvergeTol',1e-4,@isnumeric);
addParameter(p,'bpStep',0.1,@isnumeric);
addParameter(p,'bpMaxIter',500,@isnumeric);
addParameter(p,'gradStep',0.025,@isnumeric);
addParameter(p,'GsConvergeTol',1e-4,@isnumeric);

parse(p,objF,iniPt,varargin{:});

%% BFGS
interNo = 1;

% Step 1: Starting point
x0 = p.Results.iniPt;
grad0 = forwardGradient(p.Results.objF,p.Results.iniPt,p.Results.gradStep);
H0 = eye(numel(p.Results.iniPt));

% Step 2: Convergence check
while norm(grad0) >= p.Results.bfgsConvergeTol && interNo < p.Results.maxIter
    
    % Step 3: Solve set of linear equations
    direction = -H0\grad0(:);
    direction = reshape(direction,size(p.Results.iniPt,1),[]);
    % Step 4a: Create bounds for alpha_star by using bounding phase
    [alphaLeft,alphaRight] = boundingPhase(p.Results.objF,x0,direction,p.Results.bpStep);
    
    % Step 4b: Use golden section to get alpha_star
    alphaStar = goldenSection(p.Results.objF,x0,direction,alphaLeft,alphaRight,...
        p.Results.GsConvergeTol);
    
    % Step 5: Get X_new
    x1 = enforceLimits(x0 + alphaStar*direction);
    grad1 = forwardGradient(p.Results.objF,x1,p.Results.gradStep);
    
    % Step 6: Update H
    P = x1 - x0;
    Y = grad1 - grad0;
    D = (Y(:)*Y(:)')/(Y(:)'*P(:));
    E = (grad0(:)*grad0(:)')/(grad0(:)'*direction(:));
    
    H0 = H0 + D + E;
    x0 = x1;
    
    if norm(grad1-grad0) < 0.5
        H0 = eye(numel(p.Results.iniPt));
    end
    
    grad0 = grad1;
    interNo = interNo+1;
    
end

optPts = x0;
fMin = p.Results.objF(optPts);

%% secondary functions
%%%% bounding phase
    function [alphaLeft,alphaRight] = boundingPhase(objF,iniPt,direction,stepSize)
        
        maxK = p.Results.bpMaxIter;
        fVal = NaN(1,maxK);
        alpha = NaN(1,maxK);
        
        if iniPt == p.Results.lb
            fVal(1) = objF(iniPt);
            fVal(2) = objF(enforceLimits(iniPt + direction*stepSize));
            fVal(3) = objF(enforceLimits(iniPt + 2*direction*stepSize));
        elseif iniPt == p.Results.ub
            fVal(1) = objF(enforceLimits(iniPt - 2*direction*stepSize));
            fVal(2) = objF(enforceLimits(iniPt - direction*stepSize));
            fVal(3) = objF(iniPt);
        else
            fVal(1) = objF(enforceLimits(iniPt - direction*stepSize));
            fVal(2) = objF(iniPt);
            fVal(3) = objF(enforceLimits(iniPt + direction*stepSize));
        end
        
        % check direction of alha
        if fVal(1) > fVal(2) && fVal(2) > fVal(3)
            cs = 1;
        elseif fVal(1) < fVal(2) && fVal(2) < fVal(3)
            stepSize = -1*stepSize;
            cs = 1;
        else
            cs = 2;
        end
        
        % initialize alpha
        alpha(1:3) = stepSize*[-1;0;1];
        
        switch cs
            case 1
                k = 1;
                while fVal(k+2) <= fVal(k+1) && k<=maxK-2
                    k = k+1;
                    alpha(k+2) = alpha(k+1) + (2^(k-2))*stepSize;
                    fVal(k+2) = objF(enforceLimits(iniPt + direction*alpha(k+2)));
                    
                end
                alphaLeft = alpha(k);
                alphaRight = alpha(k+2);
                
                if k == maxK-2
                    error('bounding phase failed');
                end
            case 2
                alphaLeft = -stepSize;
                alphaRight = stepSize;
        end
    end

%%%% limit checker
    function X = enforceLimits(X)
        belowLim = X<p.Results.lb;
        X(belowLim) = p.Results.lb(belowLim);
        aboveLim = X>p.Results.ub;
        X(aboveLim) = p.Results.ub(aboveLim);
    end

%%%% golden section
    function alphaStar = goldenSection(objF,iniPt,direction,alphaLeft,alphaRight,...
            convergeTol)
        
        % initial length
        LStart = abs(alphaRight - alphaLeft);
        L(1) = LStart;
        
        % golden ratio
        tau = 0.381966;
        tauI = 1 - tau;
        
        alphaOne = tauI*alphaLeft + tau*alphaRight;
        alphaTwo = tau*alphaLeft + tauI*alphaRight;
        
        f1 = objF(enforceLimits(iniPt + alphaOne*direction));
        f2 = objF(enforceLimits(iniPt + alphaTwo*direction));
        
        k = 1;
        
        while L/LStart > convergeTol && k<1000
            
            if f2 > f1
                
                alphaRight = alphaTwo;
                
                alphaTwo = alphaOne;
                f2 = f1;
                
                alphaOne = tauI*alphaLeft + (tau*alphaRight);
                f1 = objF(enforceLimits(iniPt + (alphaOne*direction)));
                
                L = abs(alphaRight - alphaLeft);
                
                k = k+1;
                
            else
                
                alphaLeft = alphaOne;
                
                alphaOne = alphaTwo;
                f1 = f2;
                
                alphaTwo = (tau*alphaLeft) + tauI*alphaRight;
                f2 = objF(enforceLimits(iniPt + (alphaTwo*direction)));
                
                L = abs(alphaRight - alphaLeft);
                
                k = k+1;
                
            end
        end
        
        alphaStar = (alphaLeft + alphaRight)/2;
        
    end

%%%% forward gradient
    function grad = forwardGradient(objF,iniPt,dX)
        
        grad = NaN*iniPt;
        nEl = numel(iniPt);
        f0 = objF(iniPt);
        fNext = NaN(nEl,1);
        xNext = iniPt;
        
        for ii = 1:nEl
            xNext(ii) = xNext(ii) + dX;
            fNext(ii) = objF(enforceLimits(xNext));
            grad(ii) = (fNext(ii) - f0)/dX;
            xNext = iniPt;
        end
    end

end


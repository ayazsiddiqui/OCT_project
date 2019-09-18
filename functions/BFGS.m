function [optPts,fMin] = BFGS(objF,iniPt,lb,ub,varargin)
%BFGS Summary of this function goes here
%   Detailed explanation goes here


p = inputParser;
addRequired(p,'objF');
addRequired(p,'iniPt',@isnumeric);
addRequired(p,'lb',@isnumeric);
addRequired(p,'ub',@isnumeric);
addParameter(p,'maxIter',40,@isnumeric);
addParameter(p,'bpPerc',0.1,@isnumeric);



parse(p,objF,iniPt,lb,ub,varargin{:});

direction = [1;1];

[xLeft,xRight] = boundingPhase(p.Results.objF,p.Results.iniPt,direction,1);


%% secondary functions

%%%% bounding phase
    function [alphaLeft,alphaRight] = boundingPhase(objF,iniPt,direction,stepSize)
        
        maxK = 100;
        fVal = NaN(1,maxK);
        alpha = NaN(1,maxK);
        
        if iniPt == lb
            fVal(1) = objF(iniPt);
            fVal(2) = objF(enforceLimits(iniPt + direction*stepSize));
            fVal(3) = objF(enforceLimits(iniPt + 2*direction*stepSize));
        elseif iniPt == ub
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

% limit checker
    function X = enforceLimits(X)
        belowLim = X<lb;
        X(belowLim) = lb(belowLim);
        aboveLim = X>ub;
        X(aboveLim) = ub(aboveLim);
    end

% golden section
    function finVal = goldenSection(objF,iniPt,direction,alphaLeft,alphaRight)
        
        % initial length
        LStart = abs(alphaRight - alphaLeft);
        L(1) = LStart;
        
        % golden value
        tau = 0.381966;
        tauI = 1 - tau;
        
        k = 1;
        
        
        
    end





end


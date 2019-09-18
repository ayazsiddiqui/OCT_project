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



[xLeft,xRight] = boundingPhase(p.Results.objF,p.Results.iniPt,p.Results.bpPerc*(p.Results.ub - p.Results.lb),...
    p.Results.lb,p.Results.ub);


%% secondary functions

%%%% bounding phase
    function [xLeft,xRight] = boundingPhase(objF,iniPt,stepSize,loLim,hiLim)
        
        fMid = objF(iniPt);
        fLeft = objF(iniPt - stepSize);
        fRight = objF(iniPt + stepSize);
        
        if fLeft > fMid && fMid > fRight
            cs = 1;
        elseif fLeft < fMid && fMid < fRight
            stepSize = -1*stepSize;
            cs = 1;
        else
            cs = 2;
        end
        
        switch cs
            case 1
                X = NaN([size(iniPt) 11]);
                fVal = NaN(1,11);
                X(:,:,1) = iniPt;
                fVal(1) = fMid;
                for k = 1:10
                    X(:,:,k+1) = X(:,:,k) + (2^(k-1))*stepSize;
                    X(:,:,k+1) = enforceLimits(X(:,:,k+1),loLim,hiLim);
                    fVal(k+1) = objF(X(:,:,k+1));
                    
                    if fVal(k+1) >= fVal(k)
                        xLeft = X(:,:,k-1);
                        xRight = X(:,:,k+1);
                        break
                    end
                end
                
            case 2
                xLeft = iniPt-stepSize;
                xRight = iniPt + stepSize;
                xLeft = enforceLimits(xLeft,loLim,hiLim);
                xRight = enforceLimits(xRight,loLim,hiLim);

                
        end
        
        
    end

    function X = enforceLimits(X,lowLim,hiLim)
        belowLim = X<lowLim;
        X(belowLim) = lowLim(belowLim);
        aboveLim = X>hiLim;
        X(aboveLim) = hiLim(aboveLim);
    end




end


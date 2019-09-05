% set bounds on desgin
function val = calDesignBounds(dsgnPt,tau,designLimits)

lb = dsgnPt - tau;
ub = dsgnPt + tau;

lowLim = designLimits(:,1);
hiLim = designLimits(:,2);

belowLow = lb<lowLim;
aboveHi = ub>hiLim;

lb(belowLow) = lowLim(belowLow);
ub(aboveHi) = hiLim(aboveHi);

val = [lb(:) ub(:)];
end
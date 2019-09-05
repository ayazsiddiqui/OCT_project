% calculate log likelihood
function val = calcLogLikelihood(dsgnSet,dsgnFval,hyper)

Kmat = buildCovarianceMatrix(dsgnSet,dsgnSet,hyper(1),0,hyper(2:end));

y = dsgnFval;

val = 1*(-0.5*(y'/Kmat*y) - 0.5*log(det(Kmat)));
end
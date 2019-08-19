function op = logLogLikelihood(y,X,hyperParam)

sigma0 = hyperParam(1);
sigmaE = hyperParam(2);
theta = hyperParam(3:end);

nSamp = numel(y);
Kmat = zeros(nSamp,nSamp);
for ii = 1:nSamp
    for jj = ii:nSamp
        Kmat(ii,jj) = covFuncEval(X(:,ii),X(:,jj),sigma0,sigmaE,theta);
    end
end

Kmat = Kmat + Kmat' - eye(size(Kmat)).*diag(Kmat);


op = -1*((-0.5*y'*(Kmat\y)) - (0.5*log(norm(Kmat))));

end
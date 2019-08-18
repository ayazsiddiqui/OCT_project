function op = logLogLikelihood(y,X,hyperParam)

sigma0 = hyperParam(1);
theta = hyperParam(2:end);

nSamp = numel(y);
Kmat = zeros(nSamp,nSamp);
for ii = 1:nSamp
    for jj = ii:nSamp
        Kmat(ii,jj) = covFuncEval(X(:,ii),X(:,jj),theta,sigma0);
    end
end

Kmat = Kmat + Kmat' - eye(size(Kmat)).*diag(Kmat);


op = -1*((-0.5*y'*(Kmat\y)) - (0.5*log(norm(Kmat))));

end
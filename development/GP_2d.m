% clear
% clc
% format compact

% GP example as shown in mathemathicalmond 19.3 video GP_2d
% choose kernel (covariance function)
kernel = 2;
switch kernel
    case 1; k = @(x,y) 1*x'*y; % linear
    case 2; k = @(x,y) exp(-alpha*(x-y)'*(x-y));
end

% choose points which to sample
points = (0:0.05:1)';
[U,V] = meshgrid(points,points);
x = [U(:) V(:)]';
n = size(x,2);

% constrcut covaraince matrix
C = zeros(n,n);
for ii = 1:n
    for jj = 1:n
        C(ii,jj) = k(x(:,ii),x(:,jj));
    end
end

% sample from the GP at these points
u = randn(n,1);
[A,S,B] = svd(C);
z = A*sqrt(S)*u;

% plot
figure(2);
Z = reshape(z,sqrt(n),sqrt(n));
surf(U,V,Z)
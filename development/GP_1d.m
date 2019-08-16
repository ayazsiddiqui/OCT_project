% clear
% clc
% format compact

% GP example as shown in mathemathicalmond 19.3 video GP_1d
% choose kernel (covariance function)
alpha = 100;
beta = 5;
kernel = 6;
switch kernel
    case 1; k = @(x,y) alpha*x'*y; % linear
    case 2; k = @(x,y) alpha*min(x,y); % brownian
    case 3; k = @(x,y) exp(-alpha*(x-y)'*(x-y)); % squared exponential
    case 4; k = @(x,y) exp(-alpha*sqrt((x-y)'*(x-y))); % ornstein-Uhlenbeck
    case 5; k = @(x,y) exp(-alpha*sin(beta*pi*(x-y))^2); % periodic 
    case 6; k = @(x,y) exp(-alpha*min(abs(x-y),abs(x+y))^2); % symmetric
end

% choose points which to sample
x = (-1:0.005:1);
n = length(x);

% constrcut covaraince matrix
C = zeros(n,n);
for ii = 1:n
    for jj = 1:n
        C(ii,jj) = k(x(ii),x(jj));
    end
end

% sample from the GP at these points
u = randn(n,1);
[A,S,B] = svd(C);
z = A*sqrt(S)*u;

% plot
figure(2);
hold on
plot(x,z,'.-');
axis([-1,1,-inf,inf])
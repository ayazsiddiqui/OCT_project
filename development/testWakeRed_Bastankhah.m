clear
% clc
format compact
close all

%% attempt to recreate analytical wake redirection model from Bastankhah paper
d = 1;
yaw = 30*pi/180;
u = 1.5;

% locations
x0 = 0;
x = linspace(2*d,8*d,100);
y = linspace(-1.5*d,1.5*d,100);
z = 0.148;
zh = 0.148;

[X,Y] = meshgrid(x/d,y/d);

% wake growth rates
ky = 0.01;
kz = ky;
CT = 0.9;
TI = 0.1;
alpStar = 2.32;
betaStar = 0.154;

% wake widths
sigmaYbyD = (ky*(x-x0)/d + (cos(yaw)/sqrt(8)));
sigmaZbyD = (kz*(x-x0)/d + (1/sqrt(8)));

% normalized length of potential core
x0byD = (cos(yaw)*(1 + sqrt(1-CT)))/(sqrt(2)*(alpStar*TI + betaStar*(1-sqrt(1-CT))));

% wake skew angle
theta = ((0.3*yaw)/(cos(yaw)))*(1 - sqrt(1 - CT*cos(yaw)));

% normalized wake deflection\
t1 = theta*x0byD;
t2 = (theta/14.7)*sqrt(cos(yaw)/(ky*kz*CT));
t3 = 2.9 + 1.3*sqrt(1 - CT) - CT;
t4 = (1.6 + sqrt(CT))*(1.6*sqrt((8*sigmaYbyD.*sigmaZbyD)/cos(yaw)) - sqrt(CT));
t5 = (1.6 - sqrt(CT))*(1.6*sqrt((8*sigmaYbyD.*sigmaZbyD)/cos(yaw)) + sqrt(CT));

deltabyD = t1 + t2*t3*log(t4./t5);

% % % far wake velocity
for ii = 1:length(y)
    t1n(ii,:) = 1 - sqrt(1 - ((CT*cos(yaw))./(8*(sigmaYbyD.*sigmaZbyD))));
    t2n(ii,:) = exp(-0.5*(((y(ii)/d)-deltabyD)./sigmaYbyD).^2);
    t3n(ii,:) = exp(-0.5*(((z - zh)/d)./sigmaZbyD).^2);
    
    du(ii,:) = t1n(ii,:).*t2n(ii,:).*t3n(ii,:);
    
end

[midLine,yval] = max(du,[],1);
% for ii = 1:length(x)
%     f1 = figure(1);
%     plot(x(ii)-u*du(:,ii),y,'-k')
%     hold on
%     grid on
% end
% grid minor
% xlim([0 inf])
% xlabel('x/d')
% ylabel('y/d')

fg = figure;
contourf(X,Y,u*(1-du));
hold on
plot(x/d,y(yval)/d,'k');
grid on
colorbar
colormap jet
xlabel('x/d')
ylabel('y/d')
title(sprintf('yaw = %0.1f deg',yaw*180/pi))

saveas(fg,sprintf('yaw%0.0f.png',yaw*180/pi));


clear
clc
format compact

%% attempt to recreate analytical wake redirection model from Bastankhah paper
d = 1;
yaw = 30*pi/180;
u = 1.0;

% locations
x0 = 0;
x = 3:0.1:15;
y = -1.5*d:0.1:1.5*d;
z = -1.5*d:0.2:1.5*d;
zh = 0;

[X,Y,Z] = meshgrid(x,y,z);

% wake growth rates
ky = 0.022;
kz = 0.022;
CT = 0.9;
TI = 1;
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
% for ii = 1:
t1n = NaN(length(y),length(x),length(z));
t2n = NaN(length(y),length(x),length(z));
t3n = NaN(length(y),length(x),length(z));
du = NaN(length(y),length(x),length(z));

for jj = 1:length(z)
    for ii = 1:length(y)
        
        t1n(ii,:,jj) = 1 - sqrt(1 - ((CT*cos(yaw))./(8*(sigmaYbyD.*sigmaZbyD))));
        t2n(ii,:,jj) = exp(-0.5*((y(ii)-deltabyD)./sigmaYbyD).^2);
        t3n(ii,:,jj) = exp(-0.5*((z(jj) - zh)./sigmaZbyD).^2);
        
        du(ii,:,jj) = t1n(ii,:,jj).*t2n(ii,:,jj).*t3n(ii,:,jj);
    end
    
end
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

[xn,yn,zn,dun] = elimVals(X(:),Y(:),Z(:),du(:));

scatter3(xn,yn,zn,10,u*dun,'filled')
alpha(.15)
colorbar
xlabel('x/d')
ylabel('y/d')

figure
% space slices out without the data
sx = linspace(min(X(:)),max(X(:)),2);
sy = linspace(min(Y(:)),max(Y(:)),10);
sz = linspace(min(Z(:)),max(Z(:)),10);
% create the slices 
% in this example, there are 19 surfaces created
h = slice(X, Y, Z, u*du, sx, sy, sz);
% set properties for all 19 objects at once using the "set" function
set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
% set transparency to correlate to the data values.
alpha('color');
colormap(jet);


function [xn,yn,zn,dun] = elimVals(x,y,z,du)

xn = x;
yn = y;
zn = z;
dun = du;

fn = du<0.1;
xn(fn) = [];
yn(fn) = [];
zn(fn) = [];
dun(fn) = [];

end

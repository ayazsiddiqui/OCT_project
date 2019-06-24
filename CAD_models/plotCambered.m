clear
clc
close all
format compact

m_c = 0.02;
p = 0.4;
th = 0.12;

chord = 2;

a0 = 0.2969;
a1 = -0.126;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1036;


x1 = linspace(0,p,50)';
x2 = linspace(p,1,50)';

yc1 = (m_c/p^2)*(2*p.*x1 - x1.^2);
yc2 = (m_c/(1-p)^2)*(1 - 2*p + 2*p.*x2 - x2.^2);

dyc1 = (2*m_c/p^2)*(p - x1);
dyc2 = (2*m_c/(1-p)^2)*(p - x2);

yt1 = (th/0.2)*(a0*x1.^0.5 + a1*x1 + a2*x1.^2 + a3*x1.^3 + a4*x1.^4);
yt2 = (th/0.2)*(a0*x2.^0.5 + a1*x2 + a2*x2.^2 + a3*x2.^3 + a4*x2.^4);

xU1 = chord*(x1 - yt1.*sin(atan(dyc1)));
xU2 = chord*(x2 - yt2.*sin(atan(dyc2)));

yU1 = chord*(yc1 + yt1.*cos(atan(dyc1)));
yU2 = chord*(yc2 + yt2.*cos(atan(dyc2)));

xL1 = chord*(x1 + yt1.*sin(atan(dyc1)));
xL2 = chord*(x2 + yt2.*sin(atan(dyc2)));

yL1 = chord*(yc1 - yt1.*cos(atan(dyc1)));
yL2 = chord*(yc2 - yt2.*cos(atan(dyc2)));

%% plot
colors = 1/255*[0,0,0; 228,26,28; 55,126,184; 77,175,74];
lwd = 0.75;

plot([xU1;xU2],[yU1; yU2],'-','linewidth',lwd,'color',colors(1,:));
hold on
plot([xL1;xL2],[yL1; yL2],'-','linewidth',lwd,'color',colors(1,:));
plot(chord*[x1;x2],chord*[yc1; yc2],'-','linewidth',lwd,'color',colors(2,:));

grid on
axis equal







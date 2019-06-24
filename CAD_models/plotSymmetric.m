clear
clc
close all
format compact

th = 0.15;

chord = 2;

a0 = 0.2969;
a1 = -0.126;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1036;

x = linspace(0,1,100)';

xU = chord*x;
yU = (th/0.2)*(a0*x.^0.5 + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4);

xL = chord*x;
yL = -(th/0.2)*(a0*x.^0.5 + a1*x + a2*x.^2 + a3*x.^3 + a4*x.^4);


%% plot
colors = 1/255*[0,0,0; 228,26,28; 55,126,184; 77,175,74];
lwd = 0.75;

plot(xU,yU,'-','linewidth',lwd,'color',colors(1,:));
hold on
plot(xL,yL,'-','linewidth',lwd,'color',colors(1,:));

grid on
axis equal







clear
clc
format compact

%% attempt to recreate analytical wake redirection model
D = 1;
x = 0;
th = linspace(10*pi/180,30*pi/180,20);
Ct = 0.9;
beta = 0.125;

alp = ((cos(th).^2).*sin(th)*(Ct/2))./(1 + beta*(x./D)).^2;
alp = alp*180/pi;

%% figure
plot(x./D,alp)
plot(th*180/pi,alp)
grid on
hold on
xlabel('x/D')
ylabel('$\alpha~(deg)$')
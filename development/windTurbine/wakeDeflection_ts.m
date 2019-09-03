clear
clc
format compact

td = 0.15;

upstreamTurbPos = [0;0;0];
downstreamTurbPos = td*[7;0;0];
upstreamTurbFlow = 4;
upstreamTurbYaw = 25*pi/180;
nth = 30;
nd = 10;
ky = 0.022;
kz = 0.022;
CT = 0.9;

sim('wakeDeflection_th')
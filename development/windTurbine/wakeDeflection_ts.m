clear
clc
format compact

td = 4;

upstreamTurbPos = [0;0;0];
downstreamTurbPos = 4*[td;0;0];
upstreamTurbFlow = 4;
upstreamTurbYaw = 10*pi/180;
nth = 30;
nd = 10;
ky = 0.022;
kz = 0.022;
CT = 0.2;

sim('wakeDeflection_th')
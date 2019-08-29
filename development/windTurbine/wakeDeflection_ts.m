clear
clc
format compact


upstreamTurbPos = [0;0;0];
downstreamTurbPos = [4;0;0];
upstreamTurbFlow = 4;
upstreamTurbYaw = 10*pi/180;
nth = 30;
nd = 10;
td = 1;

sim('wakeDeflection_th')
clear
clc
format compact

td = 0.15;

upstreamTurbPos = [0;0;0];
downstreamTurbPos = td*[8;-0.8;0];
upstreamTurbFlow = 4;
upstreamTurbYaw = 10*pi/180;

Cp = 0.55;
rhoF = 1;
nth = 30;
nd = 10;
ky = 0.09;
kz = ky;
CT = 0.6;

sim('bayesian_th')
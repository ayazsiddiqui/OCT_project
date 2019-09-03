clear
clc
format compact

td = 0.15;

upstreamTurbPos = [0;0;0];
downstreamTurbPos = td*[7;0;0];
upstreamTurbFlow = 4;
upstreamTurbYaw = 30*pi/180;

Cp = 0.55;
rhoF = 1;
nth = 25;
nd = 15;
ky = 0.002;
kz = ky;
CT = 0.9;

sim('bayesian_th')
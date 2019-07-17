clear
clc

format compact

x1 = ones(3,4);
x2 = ones(3,4);

% testers
R1_o = [0;0;0];
Rn_o = [0;0;100];

N = 4;

Ri_o = [linspace(R1_o(1),Rn_o(1),N);...
    linspace(R1_o(2),Rn_o(2),N);...
    linspace(R1_o(3),Rn_o(3),N)];
Ri_o(:,1) = [];Ri_o(:,end) = [];


V1_o = [0;0;0];
Vn_o = [8;0;0];

Vi_o = [linspace(V1_o(1),Vn_o(1),N);...
    linspace(V1_o(2),Vn_o(2),N);...
    linspace(V1_o(3),Vn_o(3),N)];
Vi_o(:,1) = [];Vi_o(:,end) = [];


L = 90;
Fi = zeros(3,N);
Ti = zeros(1,N-1);
numNode = N;
massBdy = 1e3;
zetaTether = 0.05;
youngsTether = 3.8e9;
diaTether = 0.01;

[Fi1,Ti1]  = fcn(R1_o,V1_o,Rn_o,Vn_o,Ri_o,Vi_o,L,Fi,Ti,...
     numNode, massBdy, zetaTether, youngsTether, diaTether);
[Fi2,Ti2]  = fcn2(R1_o,V1_o,Rn_o,Vn_o,Ri_o,Vi_o,L,Fi,Ti,...
     numNode, massBdy, zetaTether, youngsTether, diaTether);




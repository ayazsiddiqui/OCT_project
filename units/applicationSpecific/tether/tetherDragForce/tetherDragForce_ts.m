clear
clc
format compact


Rn_o = [0 0 100]';
Vn_o = [0 0 0]';
R1_o = [0 0 0]';
V1_o = [0 0 0]';

N = 5;
Ri_o =  zeros(3,N-2);

mass = 100;
dia_t = 0.01;
rho = 1300;
rho_fluid = 1000;
Cd = 0.5;
flow = [1 0 0]';

for ii = 2:N-1
    Ri_o(:,ii-1) = (Rn_o - R1_o)*(ii-1)/(N-1);
    
end

Vi_o = zeros(size(Ri_o));


Cd*0.5*rho_fluid*norm(flow)^2*norm(Rn_o-R1_o)*(dia_t)
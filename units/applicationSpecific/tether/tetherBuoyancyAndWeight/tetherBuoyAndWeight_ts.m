clear
clc
format compact


N = 4;

rho_fluid = 1000;
rho_tether = 1010;

dia_t = 0.1;
g = 9.81;
L = 100;

vol = L*(pi/4)*dia_t^2;

W = vol*rho_tether;

open_system('tetherBuoyAndWeight_th')

sim('tetherBuoyAndWeight_th')


% hand_calc
vol = L*(pi/4)*dia_t^2;
hand_cal_ans = vol*g*(rho_fluid - rho_tether)/(N-2);

mi =  vol*(rho_tether)/(N-2);
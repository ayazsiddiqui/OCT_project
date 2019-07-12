clear
clc
format compact

base_mass = 1;

added_mass = [0 0 0;0 0 0; 0 0 10];

mass_mat = eye(3)*base_mass + added_mass;

ini_pos = [0; 0; 100];
ini_vel = [0; 0; 0];

MI = eye(3);

ini_eul = [0; 0; 0];
ini_OwB = [0; 0; 0];

FNetBdy = [0;0;-10];

MNetBdy = [0;0;0];

open_system('sixDoFDynamics_th');

sim('sixDoFDynamics_th');
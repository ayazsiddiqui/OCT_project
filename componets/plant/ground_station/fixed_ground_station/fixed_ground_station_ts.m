clear
clc
format compact

%% mask parameters
Izz = 1;
c_damp = 0.5;

class_thr(1).R1_g = [0; -1; 0];
class_thr(2).R1_g = [1; 0; 0];
class_thr(3).R1_g = [0; 1; 0];


f1 = [5;0;0];
f2 = [1;2;0];
f3 = [1;0;0];

t1 = cross(class_thr(1).R1_g,f1);
t2 = cross(class_thr(2).R1_g,f2);
t3 = cross(class_thr(3).R1_g,f3);

[t1 t2 t3]


test_forces = [f1 f2 f3];


% clear
clc
format compact

%% mask parameters
Izz = 1;
c_damp = 0.5;

class_thr(1).R1_g = [0; -1; 0];
class_thr(2).R1_g = [1; 0; 0];
class_thr(3).R1_g = [0; 1; 0];

class_thr(1).Rn_cm = [0; -0.5; -1];
class_thr(2).Rn_cm = [1; 0; -1];
class_thr(3).Rn_cm = [0; 0.5; -1];


f1 = [5;2;0];
f2 = [0;2;0];
f3 = [1;0;0];

tg1 = cross(class_thr(1).R1_g,f1);
tg2 = cross(class_thr(2).R1_g,f2);
tg3 = cross(class_thr(3).R1_g,f3);

tf1 = cross(class_thr(1).Rn_cm,f1);
tf2 = cross(class_thr(2).Rn_cm,f2);
tf3 = cross(class_thr(3).Rn_cm,f3);

FG = [tg1 tg2 tg3]


TF = [tf1 tf2 tf3]


test_forces = [f1 f2 f3];

parseLogsout

tsc.F_out;
tsc.M_out;





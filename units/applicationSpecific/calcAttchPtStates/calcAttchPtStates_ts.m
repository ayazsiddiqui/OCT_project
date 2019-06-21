clear
clc
format compact


%% tether parameters
some_class(1).tetherAttchPt = [-1; -1; 0];
some_class(2).tetherAttchPt = [1; 0; 0];
some_class(3).tetherAttchPt = [0; 1; 0];

some_class = reshape(some_class,1,[]);

% test signals
pos_test = [1;0;0];
vel_test = [0;5;0];
euler_test = [0;0;0];
angVel_test = [1;0;0];



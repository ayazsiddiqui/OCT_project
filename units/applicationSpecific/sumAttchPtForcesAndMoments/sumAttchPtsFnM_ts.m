clear
clc
format compact


%% tether parameters
some_class(1).tetherAttchPt = [-1; -1; 0];
some_class(2).tetherAttchPt = [1; 0; 0];
some_class(3).tetherAttchPt = [0; 1; 0];

some_class = reshape(some_class,1,[]);

% 
tp = plant_v2(3,2);

tp.vehicle.tetherAttchPts(1).value = [-1; -1; 0];
tp.vehicle.tetherAttchPts(2).value = [1; 0; 0];
tp.vehicle.tetherAttchPts(3).value = [0; 1; 0];

% test signals
n_t = length(tp.vehicle.tetherAttchPts);


f_test = [zeros(2,n_t); ones(1,n_t)];
euler_test = [0;0;0];

% results for euler = [0;0;0]
for ii = 1:n_t
    
    M_res(:,ii) = cross(tp.vehicle.tetherAttchPts(ii).value,f_test(:,ii));
end

M_res



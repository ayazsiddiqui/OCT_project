clear
clc
format compact


%% tether parameters
some_class(1).tetherAttchPt = [-1; -1; 0];
some_class(2).tetherAttchPt = [1; 0; 0];
some_class(3).tetherAttchPt = [0; 1; 0];

some_class = reshape(some_class,1,[]);

% test signals
n_t = length(some_class);

f_test = [zeros(2,n_t); ones(1,n_t)];
euler_test = [0;0;0];

% results for euler = [0;0;0]
for ii = 1:n_t
    
    M_res(:,ii) = cross(some_class(ii).tetherAttchPt,f_test(:,ii));
end

M_res



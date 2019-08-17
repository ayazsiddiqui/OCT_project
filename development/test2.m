clear
clc

syms a b
eq1 = a/(a+b) == 0.75;
eq2 = (a*b)/(((a+b)^2)*(a+b+1)) == (2/16);

[sol_a,sol_b] = solve([eq1,eq2],[a,b])
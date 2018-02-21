X = [1.];
df = @(t,state) dfun(t,state,@controlEQ, X);
eqC = @( parameter ) eqConstrains( parameter );

A = [ 1 ]; b = [1];
lb = [0];
ub = [(pi/2)^2];
x0 = [1];

fmincon(@objFun, x0 , [], [], [], [], [], [], @eqConstrains)
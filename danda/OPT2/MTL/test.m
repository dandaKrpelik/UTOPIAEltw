
n_cp = 10;

X = [10., zeros(1,n_cp)];
df = @(t,state) dfun(t,state,@controlEQ, X);
eqC = @( parameter ) eqConstrains( parameter );

A = [ 1 ]; b = [1];
lb = [10; -ones(n_cp,1)*(pi/2)^2];
ub = [1e12; ones(n_cp,1)*(pi/2)^2];
x0 = [1,zeros(1,n_cp)];

fmincon(@objFun, x0 , [], [], [], [], lb, ub, @eqConstrains)

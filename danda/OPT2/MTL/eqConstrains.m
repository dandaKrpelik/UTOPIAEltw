function [ C, Ceq ] = eqConstrains( X  )

Ceq = zeros(3,1);
C = 0;

alfa = X(1);
df = @(t,y) dfun(t,y, @controlEQ, X );
[t, y] = ode45(df, [0 alfa], [0 0 0 0 1]);

Ceq(1) = y(end, 2) - 0.1;
Ceq(2) = y(end, 3) - 10;
Ceq(3) = y(end,4);

end
function [u] = controlEQ(t,X)

n = length(X);

alfa = X(1);

c_vals = X(2:end);
c_points = linspace(0,1,n-1);


u = interp1(c_points, c_vals, t);
%u = sqrt(t*alfa) ;


end

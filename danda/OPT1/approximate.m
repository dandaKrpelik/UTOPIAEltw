function [ y ] = approximate( x, xs, c, phi, delta )
%APPROXIMATE Summary of this function goes here
%   Detailed explanation goes here


[n dim] = size(xs);

y=0;
for i=1:n
   r = norm(x-xs(i,:));
   y = y + c(i)*phi(r/delta); 
end


end

